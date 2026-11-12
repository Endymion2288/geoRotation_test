#!/usr/bin/env python3
"""Attach production truth to unchanged native ROOT tracks, preserving entry order."""
from array import array
from collections import Counter
import json
import math
from pathlib import Path


def accumulate(counter, prefix, reconstructed, reference):
    if not math.isfinite(reconstructed) or reconstructed == 0:
        counter[prefix+'_invalid_reco'] += 1
        return
    if not math.isfinite(reference) or reference == 0:
        counter[prefix+'_invalid_truth'] += 1
        return
    residual = abs(1000*reference/reconstructed)-1.
    counter[prefix+'_finite'] += 1
    counter[prefix+'_sum'] += residual
    counter[prefix+'_sum2'] += residual*residual
    counter[prefix+'_below_minus_half'] += residual < -.5
    counter[prefix+'_above_half'] += residual > .5
    counter[prefix+'_underflow'] += residual < -1.2
    counter[prefix+'_overflow'] += residual >= 1.2


def attach_production(directory):
    import ROOT
    directory = Path(directory)
    raw = ROOT.TFile.Open(str(directory/'track_states.root'))
    tree = raw.Get('tree')
    assert tree
    columns = ['event_nr','track_nr','majorityParticleId','t_eQOP','volume_id','layer_id',
        'chi2Sum','NDF','rMax_smt','eQOP_prt','eQOP_flt','eQOP_smt']
    tree.SetBranchStatus('*', 0)
    for column in columns:
        tree.SetBranchStatus(column, 1)
    output = ROOT.TFile(str(directory/'production.root'), 'RECREATE')
    friend = ROOT.TTree('production', 'Production truth for the same majority particle')
    qop = array('d', [math.nan])
    entry = array('q', [0])
    common = ROOT.std.vector('double')()
    link_common = ROOT.std.vector('double')()
    friend.Branch('prod_eQOP', qop, 'prod_eQOP/D')
    friend.Branch('source_entry', entry, 'source_entry/L')
    friend.Branch('prod_common', common)
    friend.Branch('link_common', link_common)
    counts = Counter()
    per_layer = {f'{s}_{l}': Counter() for s in [1,2,3] for l in [0,1,2]}
    examples = []
    maximum = 0.
    ordinal = 0
    with (directory/'production_truth.jsonl').open() as stream:
        for line in stream:
            event = json.loads(line)
            counts['source_events'] += 1
            for index, final in enumerate(event['final']):
                assert ordinal < tree.GetEntries()
                tree.GetEntry(ordinal)
                assert int(tree.event_nr) == event['event']
                assert int(tree.track_nr) == index, (ordinal, int(tree.track_nr), index)
                barcode = int(tree.majorityParticleId)
                truth = event['truth'].get(str(barcode))
                qop[0] = truth['charge']/truth['p_MeV'] if truth and truth['p_MeV'] > 0 else math.nan
                entry[0] = event['source_entry']
                common.clear()
                link_common.clear()
                selected = tree.NDF > 9 and tree.chi2Sum/tree.NDF < 20 and tree.rMax_smt < 90
                counts['writer_tracks'] += 1
                counts['github_selected_tracks'] += selected
                if selected:
                    counts['selected_majority_generator' if 0 < barcode < 200000 else 'selected_majority_secondary_or_missing'] += 1
                    counts['selected_missing_production'] += truth is None
                    counts['selected_pdg_'+str(truth['pdg']) if truth else 'selected_pdg_missing'] += 1
                for i, linked in enumerate(tree.t_eQOP):
                    linked = float(linked)
                    finite_link = math.isfinite(linked) and linked != 0
                    finite_prod = math.isfinite(qop[0]) and qop[0] != 0
                    common.push_back(qop[0] if finite_link and finite_prod else math.nan)
                    link_common.push_back(linked if finite_link and finite_prod else math.nan)
                    if not selected:
                        continue
                    station = {25:1,27:2,29:3}.get(int(tree.volume_id[i]))
                    if station is None:
                        counts['selected_other_volume_states'] += 1
                        continue
                    layer_id = int(tree.layer_id[i])
                    if layer_id not in [2,4,6]:
                        counts['selected_other_layer_states'] += 1
                        continue
                    layer = layer_id//2-1
                    m = per_layer[f'{station}_{layer}']
                    m['states'] += 1
                    m['finite_link'] += finite_link
                    m['finite_production'] += finite_prod
                    m['common_finite'] += finite_link and finite_prod
                    m['production_only'] += finite_prod and not finite_link
                    m['link_only'] += finite_link and not finite_prod
                    for component in ['prt','flt','smt']:
                        reconstructed = float(getattr(tree, 'eQOP_'+component)[i])
                        accumulate(m, component+'_native', reconstructed, linked)
                        accumulate(m, component+'_production', reconstructed, qop[0])
                        if finite_link and finite_prod:
                            accumulate(m, component+'_common_native', reconstructed, linked)
                            accumulate(m, component+'_common_production', reconstructed, qop[0])
                    if finite_link and finite_prod:
                        difference = abs(linked/qop[0]-1.)
                        maximum = max(maximum, difference)
                        m['qop_diff_gt_1e-6'] += difference > 1e-6
                        m['qop_diff_gt_1e-3'] += difference > 1e-3
                        m['qop_diff_gt_1pct'] += difference > .01
                        if difference > 1e-6 and len(examples) < 12:
                            examples.append(dict(source_entry=event['source_entry'], track=index, barcode=barcode,
                                station=station, layer=layer, link_qop_MeV_inv=linked,
                                production_qop_MeV_inv=qop[0], relative_qop_difference=difference))
                friend.Fill()
                ordinal += 1
    assert ordinal == tree.GetEntries() == friend.GetEntries()
    friend.Write()
    output.Close()
    raw.Close()
    result = dict(counts=dict(counts), per_layer=per_layer,
        max_relative_link_vs_production_qop=maximum, differing_examples=examples,
        association='file + source_entry + selected-track ordinal; production is the same majority barcode',
        units='Both reference branches in MeV^-1; native reconstructed q/p in GeV^-1')
    (directory/'truth_comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    alias = directory/'reference_3station_forward.root'
    if not alias.exists():
        alias.symlink_to('track_states.root')
    return result
