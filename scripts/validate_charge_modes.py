#!/usr/bin/env python3
"""Close actual merged tracks against native forward/backward Calypso outputs."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np


def tracks(folder):
    assert (folder/'exit_code.txt').read_text().strip() == '0'
    run = json.loads((folder/'run.json').read_text())
    assert run['success'] and run['isolate_smoothing']
    with np.load(folder/'tracks.npz') as data:
        result = {key: data[key] for key in data.files}
    return result


def signatures(table):
    result = []
    for i, entry in enumerate(table['source_entry']):
        hits = table['hit_identifiers'][table['hit_offsets'][i]:table['hit_offsets'][i+1]]
        result.append((int(entry), tuple(sorted(map(int, hits))),
                       float(table['qop_GeVinv'][i]), tuple(table['upstream_raw'][i])))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    agreement, merged = {}, {}
    for variant in ['beforeGeo', 'afterGeo']:
        for dsid in ['100116', '100117']:
            base = root/'mode_smoke_v1'/variant/dsid
            checked = 0
            for name in ['events', 'truth', 'tracks', 'segments']:
                with np.load(base/'3forward'/(name+'.npz')) as current, np.load(
                        root/'charge_smoke'/variant/dsid/'off'/(name+'.npz')) as original:
                    assert current.files == original.files
                    for field in original.files:
                        np.testing.assert_array_equal(current[field], original[field],
                            err_msg=variant+'/'+dsid+'/'+name+'/'+field)
                        checked += 1
            agreement[variant+'|'+dsid] = dict(fields=checked, all_equal=True)
            for station in ['3', '4']:
                forward = Counter(signatures(tracks(base/(station+'forward'))))
                backward = Counter(signatures(tracks(base/(station+'backward'))))
                output = Counter(signatures(tracks(base/(station+'merge'))))
                counts = Counter()
                for key, n in output.items():
                    origin = 'both' if key in forward and key in backward else (
                        'forward' if key in forward else 'backward' if key in backward else 'unresolved')
                    counts[origin] += n
                    assert n <= forward[key]+backward[key], 'Merge produced an unaccounted track'
                assert not counts['unresolved']
                merged[variant+'|'+dsid+'|'+station] = dict(
                    forward=sum(forward.values()), backward=sum(backward.values()),
                    merged=sum(output.values()), origins=dict(counts),
                    exact_hit_qop_and_upstream_parameter_closure=True)
    result = dict(default_agreement=agreement, actual_merge=merged,
                  scope='20 source events per charge/version/mode; actual native MergeTracksAlg')
    (root/'mode_smoke_v1/validation.json').write_text(json.dumps(result, indent=2)+'\n')
    (root/'mode_smoke_v1/default_agreement.json').write_text(json.dumps(agreement, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
