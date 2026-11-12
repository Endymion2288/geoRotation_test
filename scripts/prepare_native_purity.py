#!/usr/bin/env python3
"""Freeze full-file native before/after runs with isolated smoothing and purity."""
import argparse
import datetime
import difflib
import json
from pathlib import Path
import shutil


def freeze_plot_macro(root):
    """Add per-curve raw Entries to display labels; keep worker macro immutable."""
    baseline = (root/'reference/adapted/ckf_residual_perLayer.C').read_text()
    modified = baseline
    for stage, hist in [('smt', 'h1'), ('flt', 'h2'), ('prt', 'h3')]:
        old = f'Form("%s {stage}",label)'
        assert modified.count(old) == 1
        modified = modified.replace(old, f'Form("%s {stage} N=%.0f",label,{hist}->GetEntries())')
    # The original macro reuses one legend across all pads. Counts differ by
    # layer, so display a separate legend for each panel.
    legend = '''  TLegend* legF = new TLegend(0.60,0.75,0.88,0.88);
  legF->SetBorderSize(0);
  legF->SetFillStyle(0);'''
    assert modified.count(legend) == 1
    modified = modified.replace(legend, '')
    modified = modified.replace('      cCompF->cd(ilay+1);',
        '      cCompF->cd(ilay+1);\n'+legend)
    assert modified.count('if(ilay==0){') == 1
    modified = modified.replace('if(ilay==0){', '{')
    target = root/'reference/plot/ckf_residual_perLayer.C'
    target.parent.mkdir(exist_ok=True)
    target.write_text(modified)
    for name in ['CommonUtils.h', 'SetStyle.hpp']:
        shutil.copy2(root/'reference/adapted'/name, target.parent/name)
    (root/'reference/plot_labels.diff').write_text(''.join(difflib.unified_diff(
        baseline.splitlines(True), modified.splitlines(True),
        fromfile='counts/ckf_residual_perLayer.C', tofile='plot/ckf_residual_perLayer.C')))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root, reference = args.output.resolve(), args.reference.resolve()
    root.mkdir(parents=True, exist_ok=False)
    inputs = json.loads((reference/'inputs.json').read_text())
    sample, = inputs['samples']
    files = sample['files']
    actual = sorted(Path(files[0]['path']).parent.glob('*.root'))
    assert [str(p) for p in actual] == [f['path'] for f in files]
    assert all(p.stat().st_size == f['size_bytes'] for p, f in zip(actual, files))
    for relative in ['configs/reference_tracking.py', 'scripts/reference_truth_observer.py']:
        target = root/'code'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo/relative, target)
    shutil.copytree(reference/'reference', root/'reference')
    macro = root/'reference/adapted/ckf_residual_perLayer.C'
    baseline = macro.read_text()
    assert baseline.count('NDF<20&&rMax_smt<90&&NDF>9"') == 1
    modified = baseline.replace('NDF<20&&rMax_smt<90&&NDF>9"',
        'NDF<20&&rMax_smt<90&&NDF>9&&final_purity>0.7"')
    macro.write_text(modified)
    (root/'reference/purity_gate.diff').write_text(''.join(difflib.unified_diff(
        baseline.splitlines(True), modified.splitlines(True),
        fromfile='baseline/ckf_residual_perLayer.C', tofile='purity/ckf_residual_perLayer.C')))
    freeze_plot_macro(root)
    (root/'software_tools').symlink_to(reference/'software_tools', target_is_directory=True)
    cases = {}
    for case, software in [('beforeGeo', 'software_before'), ('afterGeo', 'software_after')]:
        (root/software).symlink_to(reference/software, target_is_directory=True)
        tasks = []
        for old in json.loads((reference/(case+'_tasks.json')).read_text()):
            assert old['skip'] == 0 and old['events'] == 1000
            task = dict(old, native_final_purity=True,
                runner=str(root/'code/configs/reference_tracking.py'), macro=str(macro),
                output=str(root/'shards'/case/f"f{old['file_index']:02d}_e0000"))
            tasks.append(task)
        assert len(tasks) == 20
        (root/(case+'_tasks.json')).write_text(json.dumps(tasks, indent=2)+'\n')
        smoke = dict(tasks[0], events=20, output=str(root/'smoke'/case),
            baseline_macro=str(reference/'reference/adapted/ckf_residual_perLayer.C'))
        (root/(case+'_smoke_tasks.json')).write_text(json.dumps([smoke], indent=2)+'\n')
        cases[case] = dict(software=str((reference/software).resolve()),
            material=tasks[0]['material'], initial_variance_inflation=tasks[0]['initial_variance_inflation'])
    inputs.update(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        reference=str(reference), cases=cases, truth='native t_eQOP only',
        truth_branches=['native'],
        final_purity_cut='>0.7', isolate_smoothing=True,
        scope='110051 full 20 files; 1000 entries per job; no added truth uniqueness or fiducial cuts',
        directory_and_file_sizes_rechecked=True)
    (root/'inputs.json').write_text(json.dumps(inputs, indent=2)+'\n')
    print(root)


if __name__ == '__main__':
    main()
