#!/usr/bin/env python3
"""Prepare symmetric generated inputs using the actual full Calypso G4 chain."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['smoke', 'full'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'symmetric_gun_v1_code'
    if args.phase == 'smoke':
        valid = json.loads((root/'refit_control_smoke_v1/validation.json').read_text())
        assert valid['native_58_fields_exact'] and valid['all_CKF_rows_exact']
        for folder in ['configs', 'scripts']:
            (code/folder).mkdir(parents=True, exist_ok=False)
        for name in ['momentum_simulation.py', 'slope_tracking.py']:
            shutil.copy2(repo/'configs'/name, code/'configs'/name)
        for name in ['run_slope_simulation.py', 'run_slope_task.py', 'ntuple_tracking_check.py',
                     'charge_diagnostics.py']:
            shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
        output = root/'symmetric_gun_smoke_v1'
        per_shard, shards, momenta = 20, 1, [100.]
    else:
        valid = json.loads((root/'symmetric_gun_smoke_v1/validation.json').read_text())
        assert valid['all_steps_succeeded'] and valid['actual_libraries_match']
        assert valid['identical_generated_conditions'] and valid['all_versions_same_truth']
        output = root/'symmetric_gun_20k_v1'
        per_shard, shards, momenta = 2500, 8, [100., 1000.]
    output.mkdir(exist_ok=False)
    jobs, sources = [], []
    for dsid, pdg in [('100116', -13), ('100117', 13)]:
        for momentum in momenta:
            for shard in range(shards):
                name = f'{dsid}_p{momentum:g}_s{shard:02d}'
                source = dict(code=str(code), dsid=dsid, pdg=pdg, p_GeV=momentum,
                    events=per_shard, seed=63901+100*(momentum == 1000)+shard,
                    sample_entry=per_shard*shard, file_index=shard, gun_z_F_mm=0.,
                    x_mm=0., y_mm=0., xy_half_width_mm=30., theta_half_width_rad=.0005,
                    libraries=dict(beforeGeo=str(root/'overlay_refit_before/lib'),
                                   afterGeo=str(root/'overlay_afterGeo/lib')),
                    materials=dict(beforeGeo=str(root.parent.parent/'stage02_closure/material-current-source.json'),
                        afterGeo='/afs/cern.ch/work/x/xai/public/FASER/develop/run/material-maps-alma9-rotated.json'),
                    reconstruction_modes=[dict(label='beforeGeo', variant='beforeGeo',
                        options=['--isolate-smoothing', '--charge-diagnostics']),
                        dict(label='beforeGeo_regular_both', variant='beforeGeo',
                        options=['--isolate-smoothing', '--charge-diagnostics',
                                 '--regular-refit-smoothing-chart', '--regular-refit-update-chart']),
                        dict(label='afterGeo', variant='afterGeo',
                        options=['--isolate-smoothing', '--charge-diagnostics'])],
                    output=str(output/'shards'/name))
                assert Path(source['materials']['beforeGeo']).is_file()
                path = output/(name+'_source.json')
                dump(path, source)
                jobs.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_simulation.py'),
                                                          '--task', str(path)], output=source['output']))
                sources.append(source)
    dump(output/'tasks.json', jobs)
    dump(output/'manifest.json', dict(phase=args.phase, jobs=len(jobs), source_tasks=sources,
        events_per_charge_and_momentum=per_shard*shards, seed_rule='Same particle-gun random seed for the two charges; verify actual generated positions and momenta',
        geometry='Full Calypso G4 FASERNU-04 with the default OFLCOND-FASER-05 alignment configuration at simulation, digitization and reconstruction',
        scope='Symmetric launch conditions at physical ST1 plane; true apparatus/field/material are retained. No exact mirror symmetry is assumed.',
        isolation=True, intermediate_HITS_RDO='worker scratch only'))
    print(output/'tasks.json', len(jobs), 'jobs')


if __name__ == '__main__':
    main()
