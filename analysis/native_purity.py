"""Gate unchanged native writer states with actual final-track SDO purity."""
from array import array
from collections import Counter
import json
import math
from pathlib import Path

from reference_truth import accumulate


def passes_purity(majority_hits, clusters):
    """Strict >0.7, including rejection of exactly 7/10."""
    return clusters > 0 and 10*majority_hits > 7*clusters


def seed_cuts(q, tq, theta, truth_theta):
    """problem.pdf p3 strict boundaries, q/p in the same units and theta in rad."""
    finite = all(math.isfinite(v) for v in (q, tq, theta, truth_theta)) and q != 0 and tq != 0
    return (bool(finite and abs(abs(tq/q)-1.) < .1),
            bool(finite and q*tq > 0),
            bool(finite and abs(theta-truth_theta) < 1e-4))


def read_seed_conditions(path):
    """One row per source entry, with exact seedIndex; event numbers may repeat."""
    import ROOT
    source = ROOT.TFile.Open(str(path))
    assert source and not source.IsZombie(), path
    tree = source.Get('tree')
    assert tree
    tree.SetBranchStatus('*', 0)
    names = ['seedIndex', 'majorityParticleId', 'track_QoverP', 't_QoverP',
             'track_theta', 't_theta']
    for name in ['event_number', *names]:
        tree.SetBranchStatus(name, 1)
    events, counts = [], Counter()
    seed_hist = {name: [0]*102 for name in ['momentum_all', 'momentum_good', 'theta_all', 'theta_good']}
    def fill(name, value, low, high):
        if math.isfinite(value):
            index = 0 if value < low else 101 if value >= high else 1+int((value-low)/(high-low)*100)
            seed_hist[name][index] += 1
    for row in tree:
        columns = [list(getattr(row, name)) for name in names]
        assert len(set(map(len, columns))) == 1
        seeds = {}
        barcodes = Counter()
        for index, barcode, q, tq, theta, ttheta in zip(*columns):
            cut = seed_cuts(float(q), float(tq), float(theta), float(ttheta))
            assert int(index) not in seeds
            seeds[int(index)] = dict(barcode=int(barcode), good=all(cut))
            counts['seeds'] += 1
            for label, passed in zip(['momentum', 'charge', 'theta'], cut):
                counts[label+'_pass'] += passed
            counts['momentum_charge_pass'] += cut[0] and cut[1]
            counts['all_pass'] += all(cut)
            if q and tq:
                fill('momentum_all', abs(tq/q)-1., -1.2, 1.2)
                fill('theta_all', theta-ttheta, -.003, .003)
                if all(cut):
                    fill('momentum_good', abs(tq/q)-1., -1.2, 1.2)
                    fill('theta_good', theta-ttheta, -.003, .003)
            barcodes[int(barcode)] += 1
        counts['barcode_groups_with_multiple_seeds'] += sum(n > 1 for n in barcodes.values())
        events.append(dict(event=int(row.event_number), seeds=seeds))
    source.Close()
    return events, counts, seed_hist


def attach_final_purity(directory, seed_summary=None):
    import ROOT
    import numpy as np
    directory = Path(directory)
    raw = ROOT.TFile.Open(str(directory/'track_states.root'))
    tree = raw.Get('tree')
    tree.SetBranchStatus('*', 0)
    for name in ['event_nr', 'track_nr', 'majorityParticleId', 't_eQOP', 'volume_id', 'layer_id',
                 'chi2Sum', 'NDF', 'rMax_smt', 'eQOP_prt', 'eQOP_flt', 'eQOP_smt']:
        tree.SetBranchStatus(name, 1)
    if seed_summary:
        assert tree.GetBranch('seed_nr'), 'Seed selection requires the seed-lineage writer patch'
        tree.SetBranchStatus('seed_nr', 1)
    output = ROOT.TFile(str(directory/'production.root'), 'RECREATE')
    friend = ROOT.TTree('production', 'Final output-track SDO purity; native truth unchanged')
    purity = array('d', [0.])
    entry = array('q', [0])
    seed_good = array('i', [0])
    seed_good_count = array('i', [0])
    friend.Branch('final_purity', purity, 'final_purity/D')
    friend.Branch('source_entry', entry, 'source_entry/L')
    friend.Branch('seed_good', seed_good, 'seed_good/I')
    friend.Branch('seed_good_count', seed_good_count, 'seed_good_count/I')
    seed_events, seed_counts, seed_hist = ([], Counter(), {})
    if seed_summary:
        seed_events, seed_counts, seed_hist = read_seed_conditions(seed_summary)
    output.cd()
    counts = Counter()
    per_layer = {f'{s}_{l}': Counter() for s in [1, 2, 3] for l in [0, 1, 2]}
    purity_counts = Counter()
    ordinal = 0
    event_entries, event_cumulative, event_tracks, event_symmetry = [], [], [], []
    with (directory/'production_truth.jsonl').open() as stream:
        for line in stream:
            event = json.loads(line)
            if seed_summary:
                seed_event = seed_events[counts['source_events']]
                assert seed_event['event'] == event['event']
            counts['source_events'] += 1
            event_entries.append(event['source_entry'])
            symmetry = np.zeros((9, 7), dtype=np.int64)
            selected_before = counts['selected_tracks']
            github_before = counts['github_selected_tracks']
            for index, final in enumerate(event['final']):
                assert ordinal < tree.GetEntries()
                tree.GetEntry(ordinal)
                assert int(tree.event_nr) == event['event']
                assert int(tree.track_nr) == index, (ordinal, int(tree.track_nr), index)
                n, k = final['clusters'], final['majority_hits']
                assert 0 <= k <= n
                purity[0] = k/n if n else 0.
                entry[0] = event['source_entry']
                pure = passes_purity(k, n)
                assert pure == (purity[0] > .7)
                github = tree.NDF > 9 and tree.chi2Sum/tree.NDF < 20 and tree.rMax_smt < 90
                if seed_summary:
                    seed = seed_event['seeds'][int(tree.seed_nr)]
                    seed_good_count[0] = sum(s['good'] for s in seed_event['seeds'].values())
                    seed_good[0] = int(seed['good'])
                    counts['seed_final_majority_mismatches'] += seed['barcode'] != final['majority_barcode']
                else:
                    seed_good[0] = 1
                selected = github and pure and (not seed_summary or seed_good[0])
                counts['writer_tracks'] += 1
                counts['github_selected_tracks'] += github
                counts['purity_pass_tracks'] += pure
                counts['seed_good_tracks'] += seed_good[0]
                counts['baseline_selected_tracks'] += github and pure
                counts['selected_tracks'] += selected
                if seed_summary:
                    counts['selected_seed_final_majority_mismatches'] += selected and seed['barcode'] != final['majority_barcode']
                counts['selected_native_final_majority_mismatches'] += selected and (
                    int(tree.majorityParticleId) != final['majority_barcode'])
                counts['exactly_point7_tracks'] += n > 0 and 10*k == 7*n
                purity_counts[f'{k}/{n}'] += 1
                if selected:
                    for i, reference in enumerate(tree.t_eQOP):
                        station = {25: 1, 27: 2, 29: 3}.get(int(tree.volume_id[i]))
                        layer = int(tree.layer_id[i])
                        if station is None or layer not in [2, 4, 6]:
                            continue
                        counter = per_layer[f'{station}_{layer//2-1}']
                        counter['states'] += 1
                        for stage in ['prt', 'flt', 'smt']:
                            accumulate(counter, stage+'_native', float(getattr(tree, 'eQOP_'+stage)[i]),
                                float(reference))
                        q, tq = float(tree.eQOP_smt[i]), float(reference)
                        if math.isfinite(q) and math.isfinite(tq) and q != 0 and tq != 0:
                            residual = abs(1000*tq/q)-1.
                            cell = symmetry[(station-1)*3+layer//2-1]
                            cell[0] += 1
                            for ti, threshold in enumerate([.1, .2, .5]):
                                cell[1+2*ti] += residual < -threshold
                                cell[2+2*ti] += residual > threshold
                friend.Fill()
                ordinal += 1
            event_symmetry.append(symmetry)
            event_tracks.append([len(event['final']), counts['github_selected_tracks']-github_before,
                                 counts['selected_tracks']-selected_before])
            event_cumulative.append([[[per_layer[f'{s}_{l}'][stage+'_native_'+field]
                for field in ['finite', 'below_minus_half', 'below_minus_point8']]
                for s in [1, 2, 3] for l in [0, 1, 2]] for stage in ['prt', 'flt', 'smt']])
    assert ordinal == tree.GetEntries() == friend.GetEntries()
    assert not seed_summary or len(seed_events) == counts['source_events']
    assert event_entries == list(range(event_entries[0], event_entries[0]+len(event_entries)))
    # Reading the seed summary opens and closes another ROOT file, which can
    # leave gDirectory pointing at a detached PyROOT directory.
    output.cd()
    assert friend.Write() > 0
    assert not output.TestBit(ROOT.TFile.kWriteError)
    output.Close()
    raw.Close()
    cumulative = np.asarray(event_cumulative, dtype=np.int64)
    np.savez_compressed(directory/'event_counts.npz', source_entry=event_entries,
        counts=np.diff(cumulative, axis=0, prepend=np.zeros_like(cumulative[:1])),
        tracks=np.asarray(event_tracks, dtype=np.int64),
        symmetry=np.asarray(event_symmetry, dtype=np.int64),
        symmetry_fields=['finite', 'left_point1', 'right_point1', 'left_point2', 'right_point2', 'left_half', 'right_half'],
        stages=['prt', 'flt', 'smt'],
        fields=['finite', 'below_minus_half', 'below_minus_point8'])
    result = dict(counts=dict(counts), per_layer=per_layer, purity_counts=dict(purity_counts),
        association='Input file + source entry + output-track ordinal; event_nr checked but not used as identity',
        selection='GitHub cuts AND final_purity>0.7'+(' AND seed p3 unbiased condition' if seed_summary else ''),
        truth='Unchanged native t_eQOP: SiHit-linked HepMC momentum',
        purity='max SDO-associated cluster count / all SCT clusters on final output track',
        seed_condition=('abs(abs(t_QoverP/seed_QoverP)-1)<0.1; seed_QoverP*t_QoverP>0; abs(seed_theta-t_theta)<1e-4 rad' if seed_summary else None),
        seed_summary_rows=seed_counts['seeds'], seed_events=len(seed_events),
        seed_counts=dict(seed_counts), seed_histograms=seed_hist,
        seed_association='input file + source entry + CKF trackGroup/seed_nr = seed writer seedIndex')
    (directory/'native_purity.json').write_text(json.dumps(result, indent=2)+'\n')
    (directory/'reference_3station_forward.root').symlink_to('track_states.root')
    return result
