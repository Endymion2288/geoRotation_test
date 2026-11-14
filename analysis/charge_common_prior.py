#!/usr/bin/env python3
"""Reduce the finite common-prior native-Calypso component audits."""
import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


def read_audit(path):
    rows = []
    with gzip.open(path, 'rt') as stream:
        rows.extend(json.loads(line) for line in stream)
    return rows


def truth_charge(pdg):
    if pdg in (-13, -11):
        return 1
    if pdg in (13, 11):
        return -1
    return 0


def reduce_variant(folder, dsid, variant):
    rows = read_audit(folder/'prior_audit.jsonl.gz')
    events = defaultdict(lambda: dict(truth_charge=0, truth_pdg=None, components={}))
    current = None
    for row in rows:
        stage = row['stage']
        if stage == 'event':
            current = row.get('local_event_index', row['event'])
        if current is None:
            continue
        if stage == 'truth' and row.get('signal') and row['status'] == 1 and abs(row['pdg']) == 13:
            assert events[current]['truth_pdg'] is None, (dsid, current)
            events[current]['truth_pdg'] = int(row['pdg'])
            events[current]['truth_charge'] = truth_charge(int(row['pdg']))
        if 'input_component' not in row or row['input_component'] < 0:
            continue
        component = int(row['input_component'])
        item = events[current]['components'].setdefault(component,
            dict(entry=current, component=component, found=[], selected=[]))
        if stage == 'found':
            item['found'].append(row)
        elif stage == 'selected':
            item['selected'].append(row)
    output = []
    for entry, event in sorted(events.items()):
        for component, item in sorted(event['components'].items()):
            truth = event['truth_charge']
            selected = sorted(item['selected'], key=lambda r: (-int(r['measurements']),
                                                               float(r['chi2']), int(r['candidate'])))
            first = selected[0] if selected else None
            qop = float(first['bound'][4]) if first else math.nan
            sign = truth * qop > 0 if truth and math.isfinite(qop) and qop != 0 else None
            output.append(dict(dsid=dsid, variant=variant, entry=entry, component=component,
                truth_pdg=event['truth_pdg'], truth_charge=truth, found=len(item['found']),
                selected=len(selected), first_selected_qop_GeVinv=qop,
                first_selected_sign_correct=sign,
                first_selected_measurements=int(first['measurements']) if first else None,
                first_selected_chi2=float(first['chi2']) if first else None))
    return output


def summary(rows):
    groups = defaultdict(Counter)
    for row in rows:
        g = groups[(row['dsid'], row['variant'])]
        g['components'] += 1
        g['found_components'] += row['found'] > 0
        g['selected_components'] += row['selected'] > 0
        g['no_selected_components'] += row['selected'] == 0
        g['first_selected_correct'] += row['first_selected_sign_correct'] is True
        g['first_selected_wrong'] += row['first_selected_sign_correct'] is False
        g['finite_first_qop'] += math.isfinite(row['first_selected_qop_GeVinv'])
        g['selected_tracks'] += row['selected']
    return {'|'.join(k): dict(v) for k, v in sorted(groups.items())}


def paired(rows):
    lookup = {(r['dsid'], r['entry'], r['component'], r['variant']): r for r in rows}
    out = []
    for key, before in sorted(lookup.items()):
        if key[-1] != 'beforeGeo':
            continue
        after = lookup.get((*key[:-1], 'afterGeo'))
        if after is None:
            continue
        out.append(dict(dsid=before['dsid'], entry=before['entry'], component=before['component'],
            before_found=before['found'], after_found=after['found'],
            before_selected=before['selected'], after_selected=after['selected'],
            before_qop=before['first_selected_qop_GeVinv'], after_qop=after['first_selected_qop_GeVinv'],
            before_correct=before['first_selected_sign_correct'],
            after_correct=after['first_selected_sign_correct']))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    rows = []
    for dsid in ['100044', '100116', '100117']:
        for variant in ['beforeGeo', 'afterGeo']:
            rows.extend(reduce_variant(root/dsid/variant, dsid, variant))
    output.mkdir(parents=True, exist_ok=False)
    (output/'components.json').write_text(json.dumps(rows, indent=2)+'\n')
    pairs = paired(rows)
    (output/'paired.json').write_text(json.dumps(pairs, indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps(dict(
        rows=len(rows), paired_rows=len(pairs), groups=summary(rows),
        input='Actual native Calypso component audits; one equal-weight component per finite prior split',
        selected_rule='first native selected candidate, sorted by measurements descending, chi2 ascending, candidate ID',
        sign_rule='truth charge times selected q/p > 0',
        no_component_events='Native seed rows absent from the common input are retained in the audit but do not enter component denominator'), indent=2)+'\n')
    with (output/'components.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    with (output/'paired.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=pairs[0].keys())
        writer.writeheader(); writer.writerows(pairs)
    print(json.dumps(dict(rows=len(rows), paired_rows=len(pairs), groups=summary(rows)), indent=2))


if __name__ == '__main__':
    main()
