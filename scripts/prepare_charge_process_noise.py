#!/usr/bin/env python3
"""Freeze real Calypso fixed-cluster-ID and material-noise experiments."""
import argparse
import difflib
import gzip
import json
from pathlib import Path
import shutil
import subprocess

from prepare_cluster_measurements import dump, sources


COORDINATES = ['native', 'sdo_strips_fixed', 'centroid_fixed']
PROCESSES = ['material_on', 'no_scattering', 'no_energy_loss', 'no_material']


def seed_rows(folder):
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as file:
        return [row for line in file if (row := json.loads(line))['stage'] == 'seed']


def write_inputs(root, output, original):
    key = original['dsid']+f"_{original['sample_entry']:08d}"
    folder = root/'ckf_measurement_v1/reconstruction'/original['dsid']/'native'/f"e{original['sample_entry']:08d}"
    seeds = seed_rows(folder)
    path = output/'inputs'/(key+'_seeds.txt')
    lines = []
    for row in seeds:
        values = [*row['bound'], *[v for r in row['cov'] for v in r],
                  *[v for r in row['surface'] for v in r]]
        assert len(values) == 58
        lines.append(' '.join(['0', str(row['seed']), *[format(v, '.17g') for v in values]]))
    path.write_text('\n'.join(lines)+'\n')
    previous = json.loads((root/'ckf_measurement_v1/analysis_v3/full/events.json').read_text())
    row = next(r for r in previous if r['dsid'] == original['dsid'] and
               r['sample_entry'] == original['sample_entry'])
    hits = row['native_ordered_measurement_hits']
    assert len(hits) == len(set(hits))
    with gzip.open(folder/'charge_audit.jsonl.gz', 'rt') as file:
        paths = [r for line in file if (r := json.loads(line))['stage'] == 'selected_state_trace'
                 and [step['hit'] for step in r['steps'] if step['measurement_flag']] == hits]
    assert len(paths) == 1
    selected_seed = paths[0]['seed']
    chosen = [line for line, seed in zip(lines, seeds) if seed['seed'] == selected_seed]
    assert len(chosen) == 1
    (output/'inputs'/(key+'_selected_seed.txt')).write_text(chosen[0]+'\n')
    selection = output/'inputs'/(key+'_hits.txt')
    selection.write_text(' '.join(['0', str(len(hits)), *map(str, hits)])+'\n')
    return path, selection


def build(root, output, source):
    repo = Path(__file__).resolve().parents[1]
    base = Path(json.loads((output/'manifest.json').read_text())['base_software'])
    old = Path(json.loads((base/'software.json').read_text())['source'])
    relative = Path('Tracking/Acts/FaserActsKalmanFilter/src')
    increment = ''
    for name in ['CKF2.cxx', 'CKF2.h']:
        filename = str(relative/name)
        increment += 'diff --git a/'+filename+' b/'+filename+'\n'
        increment += ''.join(difflib.unified_diff((old/filename).read_text().splitlines(True),
            (source/filename).read_text().splitlines(True), fromfile='a/'+filename, tofile='b/'+filename))
    (repo/'patches/0033_fixed_hit_full_seed_input.patch').write_text(increment)
    diff = subprocess.check_output(['git', '-C', str(source), 'diff', 'HEAD'], text=True)
    (output/'source.patch').write_text(diff)
    code = output/'code'
    for folder in ['scripts', 'configs']:
        (code/folder).mkdir(parents=True, exist_ok=False)
    for name in ['build_native_calypso.py', 'run_slope_task.py', 'run_charge_replay.py',
                 'validate_charge_audit.py', 'audit_reduce.py', 'charge_audit.py',
                 'slope_stage_extract.py', 'ntuple_tracking_check.py', 'cluster_coordinate_io.py',
                 'prepare_cluster_measurements.py', 'prepare_charge_process_noise.py']:
        shutil.copy2(repo/'scripts'/name, code/'scripts'/name)
    shutil.copy2(repo/'configs/slope_tracking.py', code/'configs/slope_tracking.py')
    dump(output/'build_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'scripts/build_native_calypso.py'), '--source', str(source),
        '--output', str(output/'software'), '--jobs', '4', '--source-patch', str(output/'source.patch')],
        output=str(output/'jobs/build'))])


def make_task(root, output, original, coordinate, process, selection=False, full=True, label=None):
    key = original['dsid']+f"_{original['sample_entry']:08d}"
    options = ['--no-seed-writer', '--actual-covariance', '--isolate-smoothing',
               '--stage-audit', '--charge-audit', '--cluster-coordinate-audit']
    if full:
        options += ['--full-seed-input-file', str(output/'inputs'/(key+('_selected_seed.txt' if selection else '_seeds.txt')))]
    if selection:
        options += ['--source-link-selection-file', str(output/'inputs'/(key+'_hits.txt')),
                    '--run-only-input-seeds']
    if coordinate == 'sdo_strips_fixed':
        options += ['--sdo-supported-strip-position']
    elif coordinate == 'centroid_fixed':
        options += ['--truth-cluster-position']
    if process in ['no_scattering', 'no_material']:
        options += ['--ckf-no-scattering']
    if process in ['no_energy_loss', 'no_material']:
        options += ['--ckf-no-energy-loss']
    mode = label or coordinate+'_'+process
    folder = output/'reconstruction'/original['dsid']/mode/f"e{original['sample_entry']:08d}"
    task = dict(original, original_output=original['output'], output=str(folder),
        runner=str(output/'code/configs/slope_tracking.py'), options=options,
        coordinate_mode=coordinate, process_mode=process, fixed_source_link_ids=selection,
        full_seed_input=full)
    if label == 'default':
        worker = 'run_charge_replay.py'
    else:
        worker = 'run_slope_task.py'
    path = output/'tasks'/(key+'_'+mode+'.json')
    dump(path, task)
    return dict(kind='analysis', command=['python3', str(output/'code/scripts'/worker),
                '--task', str(path)], output=str(folder))


def smoke(root, output):
    assert (output/'jobs/build/exit_code.txt').read_text().strip() == '0'
    overlay = output/'overlay'
    links = {}
    for folder in ['lib', 'python']:
        (overlay/folder).mkdir(parents=True, exist_ok=False)
        for library in [output/'software', root/'observer_hitcontrol_v2', root/'sdo_strip_coordinate_v1/software']:
            for item in (library/folder).iterdir():
                target = overlay/folder/item.name
                assert not target.exists()
                target.symlink_to(item.resolve(), target_is_directory=item.is_dir())
                links[str(target.relative_to(overlay))] = str(item.resolve())
    dump(overlay/'software.json', dict(files=links, actual_calypso=True))
    for folder in ['inputs', 'tasks']:
        (output/folder).mkdir()
    selected = []
    original = sources(root)
    for dsid in ['100116', '100117']:
        for wrong in [True, False]:
            selected.append(next(r for r in original if r['dsid'] == dsid and
                                 bool(r['replay_targets'][0]['wrong']) == wrong))
    tasks = []
    for row in selected:
        write_inputs(root, output, row)
        tasks += [make_task(root, output, row, 'native', 'material_on', full=False, label='default'),
                  make_task(root, output, row, 'native', 'material_on', label='identity_full'),
                  make_task(root, output, row, 'sdo_strips_fixed', 'material_on', label='identity_sdo')]
    dump(output/'smoke_tasks.json', tasks)
    dump(output/'smoke_validation_tasks.json', [dict(kind='analysis', command=['python3',
        str(output/'code/scripts/prepare_charge_process_noise.py'), '--root', str(root),
        '--phase', 'validate-smoke'], output=str(output/'jobs/validate_smoke'))])


def validate_smoke(root, output):
    from validate_charge_audit import compare_tables
    checks = []
    for job in json.loads((output/'smoke_tasks.json').read_text()):
        folder = Path(job['output'])
        assert (folder/'exit_code.txt').read_text().strip() == '0'
        task = json.loads((folder/'source_task.json').read_text())
        mode = folder.parent.name
        reference = (root/'sdo_strip_coordinate_v1/reconstruction'/task['dsid']/'sdo_strips_fixed'/folder.name
                     if mode == 'identity_sdo' else
                     root/'ckf_measurement_v1/reconstruction'/task['dsid']/'native'/folder.name)
        fields = compare_tables(folder, reference)
        assert len(fields) == 58
        left, right = seed_rows(folder), seed_rows(reference)
        assert len(left) == len(right)
        for a, b in zip(left, right):
            for key in ['seed', 'bound', 'cov', 'surface', 'position_A_mm', 'momentum_A_GeV']:
                assert a[key] == b[key], (folder, key)
        checks.append(dict(output=str(folder), reference=str(reference), fields_exact=58,
                           complete_seed_exact=True))
    dump(output/'smoke_validation.json', dict(all_exact=True, jobs=len(checks), checks=checks))


def controls(root, output):
    assert json.loads((output/'smoke_validation.json').read_text())['all_exact']
    # Six original pure wrong tracks and all 60 original normal controls.
    jobs = []
    for row in sources(root):
        write_inputs(root, output, row)
        for coordinate in COORDINATES:
            for process in PROCESSES:
                jobs.append(make_task(root, output, row, coordinate, process, selection=True))
    assert len(jobs) == 792
    dump(output/'all_tasks.json', jobs)
    smoke_keys = set()
    original = sources(root)
    for dsid in ['100116', '100117']:
        for wrong in [True, False]:
            row = next(r for r in original if r['dsid'] == dsid and
                       bool(r['replay_targets'][0]['wrong']) == wrong)
            smoke_keys.add((dsid, row['sample_entry']))
    smoke_jobs = [job for job in jobs if (Path(job['output']).parent.parent.name,
                  int(Path(job['output']).name[1:])) in smoke_keys]
    assert len(smoke_jobs) == 48
    dump(output/'fixed_smoke_tasks.json', smoke_jobs)
    done = {r['output'] for r in smoke_jobs}
    dump(output/'remaining_tasks.json', [r for r in jobs if r['output'] not in done])
    repo = Path(__file__).resolve().parents[1]
    code = output/'analysis_code_v1'
    for folder in ['scripts', 'analysis']:
        shutil.copytree(repo/folder, code/folder, ignore=shutil.ignore_patterns('__pycache__'))
    for label in ['smoke', 'full']:
        command = ['python3', str(code/'analysis/charge_process_noise.py'), '--root', str(root),
                   '--output', str(output/'analysis_v1'/label)]
        if label == 'smoke':
            command += ['--smoke']
        dump(output/(label+'_analysis_tasks.json'), [dict(kind='analysis', command=command,
            output=str(output/'jobs'/('analysis_v1_'+label)))])
    figures = output/'figures_v1'
    figures.mkdir()
    pdf_tools = '/eos/home-x/xcheng/FASER/geoRotation_test_results/stage05_geometry_comparison/software_tools'
    for phase in ['plot', 'validate']:
        dump(figures/(phase+'_tasks.json'), [dict(kind='analysis', command=['python3',
            str(code/'analysis/charge_process_noise_figures.py'), '--phase', phase,
            '--input', str(output/'analysis_v1/full'), '--output', str(figures/'results'),
            '--pdf-tools', pdf_tools], output=str(figures/'jobs'/phase))])
    dump(figures/'report_tasks.json', [dict(kind='analysis', command=['python3',
        str(code/'analysis/charge_report.py'), '--figures', str(root/'report_v9'),
        str(figures/'results'), '--output', str(root/'report_v10'), '--pdf-tools', pdf_tools],
        output=str(figures/'jobs/report'))])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--phase', choices=['build', 'smoke', 'validate-smoke', 'controls'], required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    output = root/'fixed_hit_process_v1'
    if args.phase == 'build':
        build(root, output, args.source.resolve())
    elif args.phase == 'smoke':
        smoke(root, output)
    elif args.phase == 'validate-smoke':
        validate_smoke(root, output)
    else:
        controls(root, output)
    print(output, args.phase)


if __name__ == '__main__':
    main()
