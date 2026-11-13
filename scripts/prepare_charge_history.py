#!/usr/bin/env python3
"""Freeze all available historical Calypso PHYS extraction tasks."""
import argparse
import json
from pathlib import Path
import shutil


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve();repo=Path(__file__).resolve().parents[1]
    inventory=json.loads((root/'input_inventory.json').read_text())
    code=root/'history_code';code.mkdir(exist_ok=False)
    shutil.copy2(repo/'scripts/extract_charge_history.py',code/'extract_charge_history.py')
    tasks=[]
    for label,history in inventory['historical_phys'].items():
        dsid,tag=label.split('|')
        sample=next(s for s in inventory['samples'] if s['dsid']==dsid)
        names={Path(f['path']).name:i for i,f in enumerate(sample['files'])}
        for row in history['files']:
            index=names[Path(row['path']).name.replace('PHYS.root','RDO.root')]
            source=dict(input=row['path'],events=row['entries'],dsid=dsid,tag=tag,
                file_index=index,output=str(root/'history'/dsid/tag/f'f{index:02d}'))
            path=root/f'history_source_{dsid}_{tag}_{index:02d}.json'
            path.write_text(json.dumps(source,indent=2)+'\n')
            tasks.append(dict(kind='analysis',command=['python3',str(code/'extract_charge_history.py'),
                '--task',str(path)],output=source['output']))
    (root/'history_tasks.json').write_text(json.dumps(tasks,indent=2)+'\n')
    print('Frozen historical files:',len(tasks))


if __name__=='__main__':
    main()
