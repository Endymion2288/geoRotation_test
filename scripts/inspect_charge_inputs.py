#!/usr/bin/env python3
"""Freeze all s0012 muon RDOs and the available Alma-era PHYS coverage."""
import argparse
import json
from pathlib import Path

import uproot


def inspect(path):
    with uproot.open(path) as file:
        tree_name = 'CollectionTree' if path.name.endswith('RDO.root') else 'nt'
        tree = file[tree_name]
        result = dict(path=str(path), tree=tree_name, entries=tree.num_entries,
                      size_bytes=path.stat().st_size,
                      collection_schema=tree.typenames(recursive=False))
    print(path.name, result['entries'], flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = Path('/eos/experiment/faser/data0/sim/mc24/particle_gun')
    samples, history = [], {}
    for dsid, pdg in [('100116', -13), ('100117', 13)]:
        files = [inspect(path) for path in sorted((base/dsid/'rdo/s0012').glob('*RDO.root'))]
        assert len(files) == 5, (dsid, len(files))
        samples.append(dict(dsid=dsid, pdg=pdg, production_tag='s0012',
                            files=files, total_entries=sum(f['entries'] for f in files)))
        for tag in ['s0012-dev', 's0012-r0019']:
            paths = sorted((base/dsid/'phy'/tag).glob('*PHYS.root'))
            rows = [inspect(path) for path in paths]
            present = {Path(f['path']).name.replace('PHYS.root', 'RDO.root') for f in rows}
            missing = [f['path'] for f in files if Path(f['path']).name not in present]
            history[dsid+'|'+tag] = dict(files=rows, missing_rdo_partners=missing,
                                        total_entries=sum(f['entries'] for f in rows))
    result = dict(samples=samples, historical_phys=history,
                  identity='Input production file and source-tree entry; run/event repeats are retained',
                  production_geometry='FASERNU-04', production_conditions='OFLCOND-FASER-04',
                  reconstruction_geometry='FASERNU-04', reconstruction_conditions='OFLCOND-FASER-05')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
