#!/usr/bin/env python3
"""Verify an already completed Calypso counter run after an EOS read failure."""
import argparse
import json
from pathlib import Path
import shutil

from validate_charge_audit import compare_tables


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    task = json.loads((args.source/'source_task.json').read_text())
    run = json.loads((args.source/'run.json').read_text())
    assert run['success']
    args.output.mkdir(parents=True, exist_ok=True)
    assert not (args.output/'counter_validation.json').exists()
    for label, source in [('current', args.source), ('baseline', Path(task['original_output']))]:
        folder = Path(label)
        folder.mkdir()
        for name in ['events', 'truth', 'tracks', 'segments']:
            shutil.copyfile(source/(name+'.npz'), folder/(name+'.npz'))
    fields = compare_tables(Path('current'), Path('baseline'))
    assert len(fields) == 58
    (args.output/'counter_validation.json').write_text(json.dumps(dict(all_58_fields_exact=True,
        fields=len(fields), original=task['original_output'], counter_source=str(args.source),
        original_exit_code=int((args.source/'exit_code.txt').read_text()),
        reconstruction_repeated=False, comparison='Both existing Calypso tables staged once to worker scratch'), indent=2)+'\n')
    print('Existing reconstruction verified against all 58 native fields:', args.source)


if __name__ == '__main__':
    main()
