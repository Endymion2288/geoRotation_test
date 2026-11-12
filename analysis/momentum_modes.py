#!/usr/bin/env python3
"""Summarize the real forward/backward, three/four-station and merge runs."""
import argparse
from collections import Counter,defaultdict
import gzip
import json
from pathlib import Path

import numpy as np

from pilot import event_groups,read_records,process_event,summarize
from paired_momentum import compare
from artifact_io import publish_gzip


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tasks',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    modes=['forward3','backward3','forward4','backward4','merge3','merge4']
    rows_by_mode=defaultdict(list);flows=defaultdict(Counter);extras=defaultdict(Counter)
    checks=Counter();resources=[]
    downstream_rows=defaultdict(list)
    for task in json.loads(args.tasks.read_text()):
        root=Path(task['output'])
        assert (root/'exit_code.txt').read_text().strip()=='0',root
        validation=json.loads((root/'merge_validation.json').read_text());checks.update(validation)
        resources.append(dict(variant=task['variant'],dsid=task['dsid'],resource=(root/'resource.txt').read_text()))
        for mode in modes:
            counts=Counter();rows=[];extra=extras[(task['variant'],mode,task['dsid'])]
            for ordinal,records in enumerate(event_groups(read_records(root/mode/'stages.jsonl.gz'))):
                h=records[0]
                output,flow=process_event(tuple(h[k] for k in ['run','lumi','event']),records,task['dsid'])
                counts.update(flow)
                finals={(r['seed'],r['candidate']):r for r in records if r['stage']=='final'}
                for hit in (r for r in records if r['stage']=='truth_hit'):
                    assert ((hit['sensor']>>59)&3)==hit['station'],hit
                    checks['validated_station_decodes']+=1
                for r in output:r.update(file_index=task['file_index'],source_entry=task['skip']+ordinal)
                for r in output:
                    final=r['final']
                    stations=sorted({(h>>59)&3 for h in finals[(final['seed'],final['candidate'])]['hits']}) if final else []
                    r['final_stations']=stations
                    final123=dict(final,good95=final['good95'] and {1,2,3}.issubset(stations),
                                  good100=final['good100'] and {1,2,3}.issubset(stations)) if final else None
                    downstream_rows[(task['variant'],mode)].append(dict(r,final=final123))
                rows.extend(output)
                extra['seeds_with_IFT_segment']+=sum(any(s['station']==0 for s in r.get('segments',[])) for r in records if r['stage']=='seed_circle')
                extra['final_more_than_9_layers']+=sum(r['layers']>9 for r in records if r['stage']=='final')
                for r in records:
                    if r['stage'] in ['seed','found','selected','converted','final'] and r.get('hits'):
                        pattern=','.join(map(str,sorted({(h>>59)&3 for h in r['hits']})))
                        extra[r['stage']+'_stations:'+pattern]+=1
                    if r['stage']=='final':
                        for origin in r.get('origin_collections',[]):extra['merge_origin:'+origin]+=1
            assert counts['events']==task['events']
            rows_by_mode[(task['variant'],mode)].extend(rows)
            flows[(task['variant'],mode,task['dsid'])].update(counts)
        print('Completed',task['variant'],task['dsid'],flush=True)
    result={}
    for (variant,mode),rows in rows_by_mode.items():
        result[variant+':'+mode]={}
        for dsid in sorted({r['dsid'] for r in rows}):
            selected=[r for r in rows if r['dsid']==dsid]
            z=[r['final']['position_A_mm'][0] for r in selected if r['final']]
            result[variant+':'+mode][dsid]=dict(summary=summarize(selected),
                summary_downstream123=summarize([r for r in downstream_rows[(variant,mode)] if r['dsid']==dsid]),cutflow=dict(flows[(variant,mode,dsid)]),
                mode_counts=dict(extras[(variant,mode,dsid)]),final_reference_z_F_range_mm=[min(z),float(np.median(z)),max(z)] if z else [])
    comparisons={}
    for mode in modes:
        def keyed(v):
            return {(r['dsid'],r['file_index'],r['source_entry'],r['barcode']):r for r in rows_by_mode[(v,mode)]}
        comparisons[mode]=compare(keyed('current'),keyed('strict'))
    (args.output/'summary.json').write_text(json.dumps(dict(modes=result,paired_material=comparisons,
                                                          merge_checks=dict(checks),resources=resources),indent=2)+'\n')
    local_particles=Path.cwd()/'mode_particles.jsonl.gz'
    with gzip.open(local_particles,'wt') as out:
        for (variant,mode),rows in rows_by_mode.items():
            for row in rows:out.write(json.dumps(dict(row,variant=variant,mode=mode),separators=(',',':'))+'\n')
    publish_gzip(local_particles,args.output/'particles.jsonl.gz')
    print(json.dumps(dict(merge_checks=dict(checks))))


if __name__=='__main__':main()
