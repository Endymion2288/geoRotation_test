"""Gate unchanged native writer states with actual final-track SDO purity."""
from array import array
from collections import Counter
import json
from pathlib import Path

from reference_truth import accumulate


def passes_purity(majority_hits, clusters):
    """Strict >0.7, including rejection of exactly 7/10."""
    return clusters > 0 and 10*majority_hits > 7*clusters


def attach_final_purity(directory):
    import ROOT
    import numpy as np
    directory = Path(directory)
    raw = ROOT.TFile.Open(str(directory/'track_states.root'))
    tree = raw.Get('tree')
    tree.SetBranchStatus('*', 0)
    for name in ['event_nr', 'track_nr', 'majorityParticleId', 't_eQOP', 'volume_id', 'layer_id',
                 'chi2Sum', 'NDF', 'rMax_smt', 'eQOP_prt', 'eQOP_flt', 'eQOP_smt']:
        tree.SetBranchStatus(name, 1)
    output = ROOT.TFile(str(directory/'production.root'), 'RECREATE')
    friend = ROOT.TTree('production', 'Final output-track SDO purity; native truth unchanged')
    purity = array('d', [0.])
    entry = array('q', [0])
    friend.Branch('final_purity', purity, 'final_purity/D')
    friend.Branch('source_entry', entry, 'source_entry/L')
    counts = Counter()
    per_layer = {f'{s}_{l}': Counter() for s in [1, 2, 3] for l in [0, 1, 2]}
    purity_counts = Counter()
    ordinal = 0
    event_entries, event_cumulative, event_tracks = [], [], []
    with (directory/'production_truth.jsonl').open() as stream:
        for line in stream:
            event = json.loads(line)
            counts['source_events'] += 1
            event_entries.append(event['source_entry'])
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
                selected = github and pure
                counts['writer_tracks'] += 1
                counts['github_selected_tracks'] += github
                counts['purity_pass_tracks'] += pure
                counts['selected_tracks'] += selected
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
                friend.Fill()
                ordinal += 1
            event_tracks.append([len(event['final']), counts['github_selected_tracks']-github_before,
                                 counts['selected_tracks']-selected_before])
            event_cumulative.append([[[per_layer[f'{s}_{l}'][stage+'_native_'+field]
                for field in ['finite', 'below_minus_half', 'below_minus_point8']]
                for s in [1, 2, 3] for l in [0, 1, 2]] for stage in ['prt', 'flt', 'smt']])
    assert ordinal == tree.GetEntries() == friend.GetEntries()
    assert event_entries == list(range(event_entries[0], event_entries[0]+len(event_entries)))
    friend.Write()
    output.Close()
    raw.Close()
    cumulative = np.asarray(event_cumulative, dtype=np.int64)
    np.savez_compressed(directory/'event_counts.npz', source_entry=event_entries,
        counts=np.diff(cumulative, axis=0, prepend=np.zeros_like(cumulative[:1])),
        tracks=np.asarray(event_tracks, dtype=np.int64), stages=['prt', 'flt', 'smt'],
        fields=['finite', 'below_minus_half', 'below_minus_point8'])
    result = dict(counts=dict(counts), per_layer=per_layer, purity_counts=dict(purity_counts),
        association='Input file + source entry + output-track ordinal; event_nr checked but not used as identity',
        selection='GitHub cuts AND final_purity>0.7; no unique/signal/fiducial/Good95 cut',
        truth='Unchanged native t_eQOP: SiHit-linked HepMC momentum',
        purity='max SDO-associated cluster count / all SCT clusters on final output track')
    (directory/'native_purity.json').write_text(json.dumps(result, indent=2)+'\n')
    (directory/'reference_3station_forward.root').symlink_to('track_states.root')
    return result
