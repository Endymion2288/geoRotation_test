"""Reduce full matrices on worker scratch to compact per-operation traces."""
from collections import defaultdict
import math

import numpy as np


def trace_record(records):
    first=records[0]
    result={k:first[k] for k in ['schema','run','event','lumi','collection','seed','candidate']}
    result.update(stage=first['stage']+'_trace',state=-1,trace_stage=first['stage'],steps=[])
    states=list(reversed(records))
    for i,s in enumerate(states):
        step={k:s[k] for k in ['geometry_id','state','acts_state_index','acts_previous_index','measurement_flag','outlier_flag','hole_flag']}
        step.update(z_F_mm=s['surface'][0][3],hit=s.get('hit'))
        for key in ['path_length_mm','material_full_X0','material_path_correction','material_full_thickness_mm']:
            if key in s:step[key]=s[key]
        for component in ['predicted','filtered','smoothed']:
            if component in s:
                step[component]=s[component]
                step[component+'_variance_qop']=s[component+'_cov'][4][4]
        try:
            if s['measurement_flag'] and all(k in s for k in ['predicted','filtered','measurement','projector']):
                predicted=np.array(s['predicted'],dtype=float);filtered=np.array(s['filtered'],dtype=float)
                cov=np.array(s['predicted_cov'],dtype=float);h=np.array(s['projector'],dtype=float)
                residual=np.array(s['measurement'])-h@predicted
                innovation=h@cov@h.T+np.array(s['measurement_cov']);gain=np.linalg.solve(innovation,h@cov).T
                observed=filtered[4]-predicted[4];expected=(gain@residual)[4]
                if not all(math.isfinite(v) for v in [observed,expected,*residual,*innovation.flat,*gain.flat]):raise ValueError('Nonfinite update')
                step.update(update_residual_mm=residual.tolist(),update_innovation_mm2=innovation.tolist(),
                    update_gain_qop=gain[4].tolist(),update_chi2=float(residual@np.linalg.solve(innovation,residual)),
                    update_delta_qop=float(observed),update_identity_error=float(observed-expected))
            if i+1<len(states) and 'filtered' in s and 'smoothed' in s:
                nxt=states[i+1]
                if all(k in nxt for k in ['predicted_cov','jacobian','smoothed','predicted']):
                    cf=np.array(s['filtered_cov'],dtype=float);cp=np.array(nxt['predicted_cov'],dtype=float);j=np.array(nxt['jacobian'],dtype=float)
                    gain=np.linalg.solve(cp,j@cf).T
                    delta=np.array(nxt['smoothed'])-np.array(nxt['predicted']);delta[2]=(delta[2]+math.pi)%(2*math.pi)-math.pi
                    if not np.all(np.isfinite(delta)) or not np.all(np.isfinite(gain)):raise ValueError('Nonfinite smoother')
                    step['smoothing_delta_qop']=float(s['smoothed'][4]-s['filtered'][4])
                    step['smoothing_identity_error']=float(step['smoothing_delta_qop']-(gain@delta)[4])
        except (np.linalg.LinAlgError,TypeError,ValueError):
            step['invalid_matrix_check']=True
        if i and 'filtered' in states[i-1] and 'predicted' in s and all(v is not None and math.isfinite(v) for v in [s['predicted'][4],states[i-1]['filtered'][4]]):
            step['prediction_delta_qop']=s['predicted'][4]-states[i-1]['filtered'][4]
        result['steps'].append(step)
    return result


def reduce_event(records):
    groups=defaultdict(list)
    for row in records:
        if row['stage'].endswith('_state'):
            if row['stage']!='found_state':groups[(row['stage'],row['seed'],row['candidate'])].append(row)
        else:yield row
    for group in groups.values():yield trace_record(group)
