#!/usr/bin/env python3
"""Run actual Calypso and the frozen GitHub macro, reducing ROOT on worker scratch."""
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE/'analysis' if (CODE/'analysis').exists() else CODE.parent/'analysis'))


def run(task, destination):
    started = time.monotonic()
    if 'input_entries' in task:
        import ROOT
        source = ROOT.TFile.Open(task['input'])
        assert source and not source.IsZombie()
        observed_entries = int(source.Get('CollectionTree').GetEntries())
        source.Close()
        assert observed_entries == task['input_entries']
        (destination/'input_metadata.json').write_text(json.dumps(dict(
            input=task['input'], collection_tree_entries=observed_entries,
            skip=task['skip'], events=task['events']), indent=2)+'\n')
    command = ['python3', task['runner'], '--input', task['input'], '--events', str(task['events']),
        '--skip', str(task['skip']), '--material', task['material'], '--variant', task['variant'],
        '--mask', task.get('mask', '0,1,2,3,4,5')]
    if task.get('isolate_smoothing', False):
        command.append('--isolate-smoothing')
    if task.get('writer_smoothing_check'):
        command.append('--audit-states')
    if task.get('initial_variance_inflation') is not None:
        command.extend(['--initial-variance-inflation',
            *map(str, task['initial_variance_inflation'])])
    with Path('reconstruction.log').open('w') as log:
        status = subprocess.run(['/usr/bin/time','-v','-o','resource.txt',*command],
            stdout=log, stderr=subprocess.STDOUT).returncode
    with Path('reconstruction.log').open('rb') as source, gzip.open(destination/'reconstruction.log.gz','wb') as target:
        shutil.copyfileobj(source,target)
    for name in ['run.json','resource.txt']:
        if Path(name).exists():
            shutil.copy2(name,destination/name)
    if status:
        (destination/'exit_code.txt').write_text(str(status)+'\n')
        return status
    if task.get('writer_smoothing_check'):
        from smoothing_recheck import check_writer
        check_writer(Path.cwd(), destination)
    from reference_truth import attach_production
    result = attach_production(Path.cwd())
    assert result['counts']['source_events'] == task['events']
    shutil.copy2('truth_comparison.json',destination/'truth_comparison.json')
    shutil.copy2('event_counts.npz',destination/'event_counts.npz')
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine('.L '+task['macro'])
    os.environ['REFERENCE_SCOPE_LABEL'] = f"110051 | {task['variant']} | file {task['file_index']} | {task['events']} source events"
    for name, branch in [('native','t_eQOP'),('production','prod_eQOP'),
                         ('common_native','link_common'),('common_production','prod_common')]:
        prefix = name+'_3station_counts'
        ROOT.ckf_residual_perLayer('reference_3station_forward.root',prefix,task['variant'],branch,False)
        # Keep only source-generated histograms; temporary canvases and raw trees remain on scratch.
        source = ROOT.TFile.Open(prefix+'.root')
        compact = Path(prefix+'.compact.root')
        output = ROOT.TFile(str(compact),'RECREATE')
        assert not output.IsZombie()
        for component in range(3):
            for station in [1,2,3]:
                for layer in [0,1,2]:
                    histogram = source.Get(f'h_{component}_{station}_{layer}')
                    assert histogram
                    assert histogram.Write(f'h_{component}_{station}_{layer}') > 0
        assert not output.TestBit(ROOT.TFile.kWriteError)
        output.Close()
        source.Close()
        # ROOT can report an EOS Flush failure without raising a Python exception.
        # Close on worker scratch first, then let the filesystem copy report errors.
        shutil.copy2(compact,destination/(prefix+'.root'))
        compact.unlink()
    if task.get('writer_smoothing_check'):
        from smoothing_recheck import check_smoke_histograms
        check_smoke_histograms(destination, Path(task['writer_smoothing_check']))
    resources = dict(raw_root_bytes=Path('track_states.root').stat().st_size,
        temporary_production_truth_bytes=Path('production_truth.jsonl').stat().st_size,
        worker_wall_s=time.monotonic()-started)
    (destination/'resources.json').write_text(json.dumps(resources,indent=2)+'\n')
    (destination/'task.json').write_text(json.dumps(task,indent=2)+'\n')
    (destination/'exit_code.txt').write_text('0\n')
    return 0
