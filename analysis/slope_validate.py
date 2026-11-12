#!/usr/bin/env python3
"""Check reference-figure metadata against the frozen S05 histogram summaries."""
import argparse
import json
from pathlib import Path

import numpy as np

from slope_tables import write_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--figures',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    tasks=json.loads((args.figures/'tasks.json').read_text())
    result=[]
    for task in tasks:
        command=task['command']
        source=Path(command[command.index('--input')+1])
        output=Path(command[command.index('--output')+1])
        index=json.loads((output/'index.json').read_text())
        coverage=json.loads((source/'coverage.json').read_text())
        variants={};checked=0
        for variant in ['beforeGeo','afterGeo']:
            summary=json.loads((source/(variant+'_summary.json')).read_text())
            assert coverage[variant]['complete_exact_coverage']
            assert index['source_events']==coverage[variant]['source_events']==summary['metrics']['source_events']
            for key,h in summary['histograms'].items():
                counts=np.asarray(h['counts']);edges=np.asarray(h['edges'])
                assert len(edges)==len(counts)+1 and np.all(np.diff(edges)>0)
                assert np.all(counts>=0) and np.all(counts==np.floor(counts))
                assert np.array_equal(h['sumw'],counts) and np.array_equal(h['sumw2'],counts)
                assert h['finite']+h['invalid']==h['raw_entries']
                assert int(counts.sum())+h['underflow']+h['overflow']==h['finite']
                checked+=1
            for plane in ['st1','VetoNu','sensor']:
                h=summary['histograms'][plane+'|all|physical']
                assert h['raw_entries']==summary['metrics']['long_tracks']
                assert h['invalid']==summary['metrics']['propagation'][plane]['invalid_slope']
            variants[variant]=summary
        common=set(variants['beforeGeo']['histograms'])&set(variants['afterGeo']['histograms'])
        for key in common-{'literal_github_macro'}:
            assert np.array_equal(variants['beforeGeo']['histograms'][key]['edges'],
                                  variants['afterGeo']['histograms'][key]['edges']),key
        result.append(dict(dsid=index['dsid'],source_events=index['source_events'],
            histograms_checked=checked,figure_count=index['figure_count'],
            bin_counts_flows_invalid_and_sumw2_consistent=True,common_physical_bin_edges=True,
            exact_source_coverage=True,source=str(source),figures=str(output)))
    write_json(args.output/'count_validation.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
