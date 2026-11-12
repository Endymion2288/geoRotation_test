#!/usr/bin/env python3
"""Independently check TTree scalar/vector broadcast and MeV/GeV units on a smoke."""
import argparse
import json
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def formula_divide(numerator, denominator):
    """TTreeFormula protects division by zero by returning zero."""
    numerator,denominator=np.broadcast_arrays(numerator,denominator)
    return np.divide(numerator,denominator,out=np.zeros_like(denominator,dtype=float),
        where=denominator!=0)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('directory',type=Path)
    args=parser.parse_args()
    directory=args.directory
    columns=['chi2Sum','NDF','rMax_smt','volume_id','layer_id','t_eQOP',
             'eQOP_prt','eQOP_flt','eQOP_smt']
    with uproot.open(directory/'track_states.root') as file:
        tracks=file['tree'].arrays(columns)
    with uproot.open(directory/'production.root') as file:
        truth=file['production'].arrays(['prod_eQOP','prod_common','link_common'])
    selected=(tracks.NDF>9)&(tracks.chi2Sum/tracks.NDF<20)&(tracks.rMax_smt<90)
    checks=[]
    edges=np.linspace(-1.2,1.2,101)
    for name,reference in [('native',tracks.t_eQOP),('production',truth.prod_eQOP),
                           ('common_native',truth.link_common),('common_production',truth.prod_common)]:
        with uproot.open(directory/(name+'_3station_counts.root')) as file:
            for component,stage in enumerate(['smt','flt','prt']):
                reconstructed,reference_vector=ak.broadcast_arrays(tracks['eQOP_'+stage],reference)
                for station,volume in [(1,25),(2,27),(3,29)]:
                    for layer in range(3):
                        choose=selected&(tracks.volume_id==volume)&(tracks.layer_id==2*(layer+1))
                        reco=ak.to_numpy(ak.flatten(reconstructed[choose])).astype(np.float64)
                        ref=ak.to_numpy(ak.flatten(reference_vector[choose])).astype(np.float64)
                        with np.errstate(divide='ignore',invalid='ignore'):
                            p_truth=formula_divide(1.,np.abs(ref))/1000.
                            residual=formula_divide(formula_divide(1.,np.abs(reco))-p_truth,p_truth)
                        # ROOT's axis maps NaN/+inf to overflow. This preserves the original Draw.
                        bins=np.searchsorted(edges,residual,side='right')
                        counted=np.bincount(bins,minlength=102)
                        hist=file[f'h_{component}_{station}_{layer}']
                        observed=hist.values(flow=True)
                        assert np.array_equal(counted,observed),(name,stage,station,layer,counted-observed)
                        checks.append(dict(truth=name,stage=stage,station=station,layer=layer,
                            states=len(residual),finite_residual=int(np.isfinite(residual).sum()),
                            flow_bins_exactly_equal=True))
    output=dict(checks=checks,histograms_checked=len(checks),
        purpose='Independent smoke check of ROOT alias broadcasting, units and NaN flow bins')
    (directory/'independent_histogram_validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(dict(histograms_checked=len(checks),all_flow_bins_equal=True)))


if __name__=='__main__':
    main()
