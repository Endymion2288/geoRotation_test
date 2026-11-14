#!/usr/bin/env python3
"""Pair actual ACTS backward operations after exact observer replay gates."""
import argparse
from collections import Counter, defaultdict
import gzip
from itertools import zip_longest
import json
from pathlib import Path

import numpy as np

from charge_component_refit import audit, components, first, metrics
from charge_refit_input import observer_gate, refits
from slope_tables import write_json
from validate_charge_components import validate


D = np.diag([1., 1., 1., -1., 1., 1.])


def gate(root, output, mode):
    checks = []
    for point in json.loads((root/'manifest.json').read_text())['points']:
        old, new = Path(point['reference']), Path(point[mode])
        check = observer_gate(components(old), components(new))
        assert refits(old) == refits(new)
        # Check the complete saved stream as well, including truth, found,
        # component errors/end markers and all refit records.
        count = 0
        for a, b in zip_longest(audit(old/'prior_audit.jsonl.gz'), audit(new/'prior_audit.jsonl.gz')):
            assert a == b, (point['dsid'], point['variant'], mode, count)
            count += 1
        array_fields=0
        archives=['events','truth','tracks','segments']
        if (old/'refit_outcomes.npz').exists():
            archives.append('refit_outcomes')
        for name in archives:
            with np.load(old/(name+'.npz'),allow_pickle=False) as a, np.load(new/(name+'.npz'),allow_pickle=False) as b:
                assert a.files==b.files
                for field in a.files:
                    np.testing.assert_array_equal(a[field],b[field],err_msg=name+':'+field)
                    array_fields+=1
        checks.append(dict(dsid=point['dsid'], variant=point['variant'], mode=mode,
                           **check, complete_audit_rows=count, runtime=validate(new),
                           refit_records_exact=True, complete_audit_stream_exact=True,
                           complete_observer_arrays_exact=True,observer_array_fields=array_fields))
        if mode == 'off':
            assert not (new/'smoothing_operations.jsonl.gz').exists()
        else:
            assert (new/'smoothing_operations.jsonl.gz').is_file()
    output.mkdir(parents=True, exist_ok=False)
    write_json(output/'validation.json', dict(mode=mode, all_exact=True, checks=checks,
        scope='Actual compiled Calypso full-stream replay against the frozen analytic-input reference.'))
    print(json.dumps(dict(mode=mode, exact_replays=len(checks))))


def common(row):
    result = {}
    for key, value in row.items():
        # Native output snapshots are in the original before-F/after-A
        # charts; only the calculate() operand/result buffers share the
        # regular-before versus native-after charts compared here.
        if isinstance(value, list) and not key.startswith('native_'):
            array = np.asarray(value, dtype=float)
            if not row['regular']:
                if array.shape == (6,):
                    array = D @ array
                    array[3] += np.pi/2.
                elif array.shape == (6,6):
                    array = D @ array @ D
            result[key] = array
    return result


def scaled(a, b, diagonal_a=None, diagonal_b=None):
    da = np.diag(a) if diagonal_a is None else diagonal_a
    db = np.diag(b) if diagonal_b is None else diagonal_b
    diagonal = np.maximum((np.abs(da)+np.abs(db))/2., 1e-30)
    return np.abs(a-b)/np.sqrt(np.outer(diagonal, diagonal))


def condition(matrix):
    diagonal = np.sqrt(np.maximum(np.abs(np.diag(matrix)), 1e-300))
    normalized = matrix/np.outer(diagonal,diagonal)
    return dict(raw=float(np.linalg.cond(matrix)), correlation_scaled=float(np.linalg.cond(normalized)))


def operation_rows(folder):
    rows = defaultdict(list)
    for row in audit(folder/'smoothing_operations.jsonl.gz'):
        assert row['raw_covariance_captured']
        assert all(np.all(np.isfinite(value)) for value in common(row).values())
        rows[row['event'],row['component'],row['candidate']].append(row)
    for path in rows.values():
        assert [r['operation'] for r in path] == list(range(len(path)))
        for a,b in zip(path,path[1:]):
            assert b['next_state_index'] == a['state_index']
            # Each next-S operand must be the exact native result from the
            # preceding backward operation (up to the existing chart map).
            np.testing.assert_allclose(common(b)['next_smoothed_mean'],
                common(a)['calculated_smoothed_mean'], rtol=2e-12, atol=2e-12)
    return rows


def analyze(root, output):
    for mode in ('off','on'):
        assert json.loads((root/('gate_'+mode+'_v1')/'validation.json').read_text())['all_exact']
    points = json.loads((root/'manifest.json').read_text())['points']
    by_sample = {p['dsid']: {} for p in points}
    for point in points:
        folder = Path(point['on'])
        by_sample[point['dsid']][point['variant']] = dict(
            rows=operation_rows(folder), tracks=components(folder), folder=str(folder))
    output.mkdir(parents=True, exist_ok=False)
    paired, summaries, paths, extremes, checks = [], {}, [], [], []
    for dsid, variants in by_sample.items():
        before, after = (variants[v] for v in ('beforeGeo','afterGeo'))
        counts = Counter()
        sample = []
        for key,a in before['tracks'].items():
            b=after['tracks'][key]
            ia, ib=first(a), first(b)
            if ia is None or ib is None:
                continue
            counts['common_selected_components'] += 1
            if a['selected'][ia]['hits'] != b['selected'][ib]['hits']:
                counts['different_hit_paths'] += 1
                continue
            pa,pb=before['rows'].get((*key,ia),[]),after['rows'].get((*key,ib),[])
            if pa:
                counts['paired_paths_with_operations']+=1
            else:
                assert len(a['states'][ia])==len(b['states'][ib])==1
                counts['terminal_only_paths']+=1
            assert [(r['geometry_id'],r['hit'],r['measurement_flag']) for r in pa] == [
                (r['geometry_id'],r['hit'],r['measurement_flag']) for r in pb]
            path=[]
            for ra,rb in zip(pa,pb):
                ca,cb=common(ra),common(rb)
                counts['paired_backward_operations'] += 1
                item=dict(dsid=dsid, entry=key[0], component=key[1],
                    candidate_before=ia, candidate_after=ib, operation=ra['operation'],
                    hit=ra['hit'], geometry_id=ra['geometry_id'],
                    measurement_flag=ra['measurement_flag'], path_length_mm=ra['path_length_mm'])
                for field in ('filtered_C','next_predicted_C','next_smoothed_C',
                              'raw_smoothed_C','validated_smoothed_C'):
                    delta=scaled(ca[field],cb[field])
                    item[field+'_difference']=dict(all=float(delta.max()), direction=float(delta[2:4,2:4].max()))
                for field in ('filtered_mean','next_predicted_mean','next_smoothed_mean','calculated_smoothed_mean'):
                    delta=np.abs(ca[field]-cb[field])
                    item[field+'_difference']=dict(direction_rad=float(delta[2:4].max()),
                        qop_GeVinv=float(delta[4]), local_mm=float(delta[:2].max()))
                item['transport_J_max_difference']=float(np.abs(ca['transport_J']-cb['transport_J']).max())
                item['gain_max_difference']=float(np.abs(ca['reconstructed_gain']-cb['reconstructed_gain']).max())
                for name, raw, c in [('before',ra,ca),('after',rb,cb)]:
                    correction=scaled(c['raw_smoothed_C'],c['validated_smoothed_C'])
                    item[name]=dict(raw_psd=raw['raw_psd'],validated_psd=raw['validated_psd'],
                        psd_correction_exactly_zero=bool(np.array_equal(c['raw_smoothed_C'],c['validated_smoothed_C'])),
                        psd_correction_scaled_max=float(correction.max()),
                        predicted_C_condition=condition(c['next_predicted_C']),
                        gain_infinity_norm=float(np.linalg.norm(c['reconstructed_gain'],ord=np.inf)))
                    counts[name+'_raw_not_psd'] += not raw['raw_psd']
                    counts[name+'_validated_not_psd'] += not raw['validated_psd']
                    counts[name+'_corrected_C'] += not item[name]['psd_correction_exactly_zero']
                sample.append(item)
                paired.append(item)
                path.append(item)
            # Descriptive, not selection thresholds. Execution order is the
            # actual backward order, so next-S is inherited from the prior op.
            for threshold in (1e-8,1e-6,1e-4):
                found=next((r for r in path if r['validated_smoothed_C_difference']['direction']>threshold),None)
                if found:
                    paths.append(dict(threshold=threshold,**found))
        summaries[dsid]=dict(counts=counts,
            fields={field:metrics([r[field]['direction'] for r in sample]) for field in
                    ('filtered_C_difference','next_predicted_C_difference','next_smoothed_C_difference',
                     'raw_smoothed_C_difference','validated_smoothed_C_difference')},
            condition={v:metrics([r[v]['predicted_C_condition']['correlation_scaled'] for r in sample])
                       for v in ('before','after')})
        if sample:
            for field in ('raw_smoothed_C_difference','validated_smoothed_C_difference'):
                worst=max(sample,key=lambda r:r[field]['direction'])
                ka=(worst['entry'],worst['component'],worst['candidate_before'])
                kb=(worst['entry'],worst['component'],worst['candidate_after'])
                extremes.append(dict(field=field, metrics=worst,
                    before_operations=before['rows'][ka],after_operations=after['rows'][kb]))
        checks.append(dict(dsid=dsid,all_logged_raw_C_captured=True,
            execution_order_and_geometry_hit_flags_match=True,source_before=before['folder'],source_after=after['folder']))
    for name,rows in [('paired_operations',paired),('first_backward_growth',paths)]:
        with gzip.open(output/(name+'.jsonl.gz'),'wt',compresslevel=4) as stream:
            for row in rows:
                stream.write(json.dumps(row,separators=(',',':'),allow_nan=False)+'\n')
    write_json(output/'extreme_operation_paths.json',extremes)
    write_json(output/'summary.json',dict(groups=summaries,
        scope='Actual Calypso ACTS smoothing operations, frozen n32 C/256 analytic input; '
              'common first-ranked selected paths. Reconstructed gain is diagnostic. '
              'Linear native-after polar-to-latitude re-expression only; no new fit.'))
    write_json(output/'validation.json',dict(checks=checks,default_and_enabled_full_stream_gates=True,
        all_raw_C_accesses_captured=True,paired_operations=len(paired)))
    print(json.dumps(dict(samples=len(summaries),paired_operations=len(paired))))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--phase',choices=['gate-off','gate-on','analyze'],required=True)
    args=parser.parse_args()
    if args.phase.startswith('gate-'):
        gate(args.input,args.output,args.phase[5:])
    else:
        analyze(args.input,args.output)


if __name__=='__main__':
    main()
