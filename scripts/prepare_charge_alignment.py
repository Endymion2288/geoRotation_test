#!/usr/bin/env python3
"""Freeze native/default and matched nominal-alignment Calypso controls."""
import argparse
import json
from pathlib import Path
import shutil


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def freeze_code(repo, code):
    for folder in ['configs', 'scripts']:
        (code/folder).mkdir(parents=True, exist_ok=False)
    for name in ['slope_tracking.py', 'momentum_simulation.py']:
        shutil.copy2(repo/'configs'/name, code/'configs'/name)
    for name in ['run_slope_task.py', 'run_slope_simulation.py', 'ntuple_tracking_check.py',
                 'charge_diagnostics.py', 'calypso_alignment.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--phase', choices=['compatibility', 'gun-smoke', 'gun-recovery'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    repo = Path(__file__).resolve().parents[1]
    code = root/'alignment_control_v1_code'
    software = root/'software_nominal_alignment_v1'
    assert (root/'nominal_alignment_build_v1_job/exit_code.txt').read_text().strip() == '0'
    libraries = {variant: str(root/('alignment_overlay_'+variant+'_v1')/'lib')
                 for variant in ['beforeGeo', 'afterGeo']}
    if args.phase == 'compatibility':
        for variant, base in [('beforeGeo', root/'overlay_refit_before'), ('afterGeo', root/'overlay_afterGeo')]:
            overlay = Path(libraries[variant]).parent
            overlay.mkdir(exist_ok=False)
            files = {}
            for folder in ['lib', 'python']:
                (overlay/folder).mkdir()
                for source in [base, software]:
                    for item in (source/folder).iterdir():
                        target = overlay/folder/item.name
                        assert not target.exists()
                        target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                        files[str(target.relative_to(overlay))] = str(item.resolve())
            dump(overlay/'software.json', dict(sources=[str(base), str(software)], files=files,
                before_after_conditions_package_original_sources_identical=True))
        freeze_code(repo, code)
        output = root/'alignment_compatibility_v1'
        output.mkdir(exist_ok=False)
        tasks, sources = [], []
        for variant in ['beforeGeo', 'afterGeo']:
            jobs = []
            for dsid in ['100044', '100116', '100117']:
                original = root/'charge_smoke'/variant/dsid/'off'
                source = json.loads((original/'source_task.json').read_text())
                source.update(runner=str(code/'configs/slope_tracking.py'), original_output=str(original),
                    output=str(output/variant/dsid),
                    options=source['options']+['--alignment-audit'])
                path = output/(variant+'_'+dsid+'.json')
                dump(path, source)
                jobs.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_task.py'),
                    '--task', str(path)], output=source['output']))
                sources.append(source)
            dump(output/(variant+'_tasks.json'), jobs)
            tasks.extend(jobs)
        dump(output/'manifest.json', dict(jobs=len(tasks), source_tasks=sources, libraries=libraries,
            control='Conditions property defaults unchanged; read-only actual/default sensor audit enabled; compare all native tables with original 20-event smoke'))
    elif args.phase == 'gun-smoke':
        valid = json.loads((root/'alignment_compatibility_v1/validation.json').read_text())
        assert valid['all_58_fields_exact'] and valid['actual_libraries_match']
        code = root/'alignment_gun_v1_code'
        freeze_code(repo, code)
        output = root/'alignment_gun_smoke_v1'
        output.mkdir(exist_ok=False)
        tasks, sources = [], []
        originals = json.loads((root/'symmetric_gun_smoke_v1/manifest.json').read_text())['source_tasks']
        for original in originals:
            for nominal in [False, True]:
                name = original['dsid']+('_nominal' if nominal else '_default')
                source = dict(original, code=str(code), libraries=libraries,
                    nominal_alignment=nominal, original_output=original['output'],
                    output=str(output/name))
                source['reconstruction_modes'] = [dict(mode, options=mode['options']+['--alignment-audit'])
                                                  for mode in original['reconstruction_modes']]
                path = output/(name+'.json')
                dump(path, source)
                tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_simulation.py'),
                    '--task', str(path)], output=source['output']))
                sources.append(source)
        dump(output/'tasks.json', tasks)
        dump(output/'manifest.json', dict(jobs=len(tasks), source_tasks=sources, libraries=libraries,
            control='Same exact launch conditions and seeds; default versus nominal alignment at full G4, digitization and all native reconstruction modes; intermediate pools remain on scratch'))
    else:
        originals = json.loads((root/'alignment_gun_smoke_v1/manifest.json').read_text())['source_tasks']
        originals = [row for row in originals if row['nominal_alignment']]
        code = root/'alignment_gun_v2_code'
        freeze_code(repo, code)
        output = root/'alignment_gun_recovery_v2'
        output.mkdir(exist_ok=False)
        tasks, sources = [], []
        for original in originals:
            previous = Path(original['output'])
            assert (previous/'exit_code.txt').read_text().strip() == '1'
            source = dict(original, code=str(code), output=str(output/previous.name), failed_output=str(previous))
            path = output/(previous.name+'.json')
            dump(path, source)
            tasks.append(dict(kind='analysis', command=['python3', str(code/'scripts/run_slope_simulation.py'),
                '--task', str(path)], output=source['output']))
            sources.append(source)
        dump(output/'tasks.json', tasks)
        dump(output/'manifest.json', dict(jobs=len(tasks), source_tasks=sources,
            recovery='Two terminal nominal jobs: G4 completed, but the audit ran after application teardown and could not retrieve GeoModel. Move the read-only audit into the first event; default jobs continue unchanged.'))
    print(output, len(tasks), 'jobs')


if __name__ == '__main__':
    main()
