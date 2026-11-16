"""Reduce actual latest-Calypso states; separate transport, update and smoothing in q/p."""
from collections import defaultdict
import gzip
import json
import math
import subprocess
from pathlib import Path

import numpy as np


SEED_FIELDS = ['seedIndex', 'majorityParticleId', 'track_loc0', 'track_loc1',
    'track_phi', 'track_theta', 'track_QoverP', 'track_time',
    't_loc0', 't_loc1', 't_phi', 't_theta', 't_QoverP', 't_time']
WRITER_SCALARS = ['event_nr', 'track_nr', 'majorityParticleId', 'chi2Sum', 'NDF', 'rMax_smt']
WRITER_VECTORS = ['volume_id', 'layer_id', 't_eQOP', 'eQOP_prt', 'eQOP_flt', 'eQOP_smt',
                  'err_eQOP_prt', 'err_eQOP_flt', 'err_eQOP_smt']
STATE_FIELDS = ['entry', 'candidate', 'seed', 'station', 'layer', 'measurement', 'z_F_mm',
    'p_production_GeV', 'native_truth_qop_GeV_inv', 'q_predicted', 'q_filtered', 'q_smoothed',
    'var_predicted', 'var_filtered', 'var_smoothed', 'previous_q_filtered',
    'dq_transport', 'dq_update', 'dq_smoothing', 'measurement_residual_mm', 'innovation_mm2',
    'gain_qop_per_mm', 'update_identity_error', 'seed_qop', 'seed_truth_qop',
    'github', 'purity', 'good_seed', 'same_seed_truth', 'material_X0', 'sigma_over_truth',
    'truth_charge', 'seed_barcode', 'track_barcode', 'predicted_corr_min_eigenvalue',
    'filtered_corr_min_eigenvalue', 'corr_qop_phi', 'variance_update_identity_error']
TRACK_FIELDS = ['entry', 'candidate', 'seed', 'seed_barcode', 'track_barcode', 'final_barcode',
    'github', 'purity', 'good_seed', 'p_production_GeV', 'seed_qop', 'seed_truth_qop']


def grouped_audit(path):
    if not path.exists():
        return
    rows = []
    with path.open() as stream:
        for line in stream:
            row = json.loads(line)
            if row['stage'] == 'event' and rows:
                yield rows
                rows = []
            rows.append(row)
    if rows:
        yield rows


def seed_good(values):
    q, tq, theta, truth_theta = values[6], values[12], values[5], values[11]
    return all(map(math.isfinite, [q,tq,theta,truth_theta])) and q*tq>0 and abs(abs(tq/q)-1)<.1 and abs(theta-truth_theta)<1e-4


def covariance_health(matrix):
    matrix = np.asarray(matrix, dtype=float)
    diag = np.diag(matrix)
    if not np.all(np.isfinite(matrix)) or np.any(diag<=0):
        return math.nan
    corr = matrix/np.sqrt(np.outer(diag,diag))
    return float(np.linalg.eigvalsh((corr+corr.T)/2)[0])


def reduce_run(directory, destination, task):
    import uproot
    directory, destination = Path(directory), Path(destination)
    with uproot.open(directory/'reference_3station_forward_seed_summary_circleFitTrackSeedTool.root') as source:
        seed_columns = source['tree'].arrays(['event_number',*SEED_FIELDS],library='np')
    seed_events, seed_rows = [], []
    for entry in range(len(seed_columns['event_number'])):
        columns=[seed_columns[name][entry] for name in SEED_FIELDS]
        assert len(set(map(len,columns)))==1
        seeds={}
        for values in zip(*columns):
            values=list(map(float,values))
            seeds[int(values[0])]=values
            seed_rows.append([task['skip']+entry,*values])
        seed_events.append((int(seed_columns['event_number'][entry]),seeds))
    assert len(seed_events)==task['events']
    with uproot.open(directory/'track_states.root') as source:
        native=source['tree'].arrays(WRITER_SCALARS+WRITER_VECTORS,library='np')
    native_entries=len(native['event_nr'])
    audits = iter(grouped_audit(directory/'momentum_audit.jsonl'))
    observed = (directory/'momentum_audit.jsonl').exists()
    writer_scalars, writer_vectors, offsets = [], [], [0]
    states_out, geometry_ids, hit_ids, tracks_out = [], [], [], []
    purity_values, seed_values_out = [], []
    ordinal = 0
    with (directory/'production_truth.jsonl').open() as stream:
        for local_entry, line in enumerate(stream):
            final_event = json.loads(line)
            event_number, seeds = seed_events[local_entry]
            assert final_event['event'] == event_number
            audit = next(audits) if observed else []
            truth = {r['barcode']:r for r in audit if r['stage']=='truth'}
            initial = {r['seed']:r for r in audit if r['stage']=='seed'}
            selected = {r['candidate']:r for r in audit if r['stage']=='selected'}
            state_map = defaultdict(list)
            for row in audit:
                assert row['source_entry'] == local_entry
                if row['stage']=='selected_state':
                    state_map[row['candidate']].append(row)
            assert not observed or len(selected) == len(final_event['final'])
            for candidate, final in enumerate(final_event['final']):
                assert ordinal<native_entries
                scalar = [float(native[name][ordinal]) for name in WRITER_SCALARS]
                assert int(scalar[0]) == event_number and int(scalar[1]) == candidate
                arrays = [np.asarray(native[name][ordinal], dtype=np.float64) for name in WRITER_VECTORS]
                assert len(set(map(len,arrays))) == 1
                writer_scalars.append(scalar)
                writer_vectors.extend(np.array(arrays).T.tolist())
                offsets.append(len(writer_vectors))
                github = scalar[4]>9 and scalar[3]/scalar[4]<20 and scalar[5]<90
                purity = final['majority_hits']/final['clusters'] if final['clusters'] else 0.
                seed = selected[candidate]['seed'] if observed else -1
                seed_values = seeds[seed] if observed else None
                good = seed_good(seed_values) if observed else False
                purity_values.append(purity);seed_values_out.append(int(good))
                if observed:
                    barcode, seed_barcode = int(scalar[2]), int(seed_values[1])
                    tr = truth.get(barcode)
                    p = np.linalg.norm(tr['momentum_F_MeV'])*.001 if tr else math.nan
                    charge = (-1. if tr['pdg']>0 else 1.) if tr else math.nan
                    tracks_out.append([task['skip']+local_entry,candidate,seed,seed_barcode,barcode,
                        final['majority_barcode'],github,purity,good,p,initial[seed]['bound'][4],seed_values[12]])
                    # Audit and ROOT writer use the same reverse-ordered measurement states.
                    reverse_states = state_map[candidate]
                    measurements = [s for s in reverse_states if s['measurement_flag']]
                    assert len(measurements) == len(arrays[0])
                    native_truth = {}
                    for j,s in enumerate(measurements):
                        assert (s['geometry_id']>>56)==int(arrays[0][j])
                        for k,name in enumerate(('predicted','filtered','smoothed')):
                            np.testing.assert_allclose(s[name][4],arrays[3+k][j],rtol=1e-6,atol=1e-12)
                        native_truth[s['state']] = arrays[2][j]*1000.
                    previous = initial[seed]['bound'][4]
                    for s in reversed(reverse_states):
                        # Hole/material states have no allocated filtered state.
                        # Their continuing mean is the predicted mean, before
                        # that surface's material interaction.
                        effective = dict(s)
                        if 'filtered' not in effective:
                            effective['filtered'] = s['predicted']
                            effective['filtered_cov'] = s['predicted_cov']
                        values = [effective.get(n,[math.nan]*6)[4] for n in ('predicted','filtered','smoothed')]
                        variances = [effective.get(n,[[math.nan]*6]*6)[4][4] for n in
                            ('predicted_cov','filtered_cov','smoothed_cov')]
                        qp,qf,qs = values
                        residual, innovation, gain = math.nan,math.nan,math.nan
                        if s.get('measurement_flag') and 'measurement' in s:
                            H = np.asarray(s['projector'],dtype=float)
                            C = np.asarray(s['predicted_cov'],dtype=float)
                            innovation = float((H@C@H.T+np.asarray(s['measurement_cov']))[0,0])
                            residual = float((np.asarray(s['measurement'])-H@np.asarray(s['predicted']))[0])
                            gain = float((C@H.T)[4,0]/innovation)
                        error = (qf-qp)-gain*residual
                        covariance_error = math.nan
                        if math.isfinite(gain):
                            covariance_error = variances[1]-(variances[0]-gain*gain*innovation)
                        cp = np.asarray(effective['predicted_cov'],dtype=float)
                        cqphi = cp[4,2]/math.sqrt(cp[4,4]*cp[2,2]) if cp[4,4]*cp[2,2]>0 else math.nan
                        geo = s['geometry_id']
                        station = {25:1,27:2,29:3}.get(geo>>56,0)
                        layer = ((geo>>36)&0xfff)//2-1
                        states_out.append([task['skip']+local_entry,candidate,seed,station,layer,
                            s['measurement_flag'],s['surface'][0][3],p,native_truth.get(s['state'],math.nan),
                            *values,*variances,previous,qp-previous,qf-qp,qs-qf,residual,innovation,gain,error,
                            initial[seed]['bound'][4],seed_values[12],github,purity,good,seed_barcode==barcode,
                            s.get('material_full_X0',0.), math.sqrt(variances[0])*p if variances[0]>=0 else math.nan,
                            charge,seed_barcode,barcode,
                            covariance_health(effective['predicted_cov']),
                            covariance_health(effective['filtered_cov']),cqphi,covariance_error])
                        geometry_ids.append(geo)
                        hit_ids.append(s.get('hit',0))
                        previous = qf
                ordinal += 1
    assert ordinal==native_entries
    assert not observed or next(audits,None) is None
    with uproot.recreate(directory/'production.root') as friend:
        friend.mktree('production',{'final_purity':'float64','seed_good':'int32'})
        if ordinal:
            friend['production'].extend(dict(final_purity=np.asarray(purity_values),
                                            seed_good=np.asarray(seed_values_out,dtype=np.int32)))
    np.savez_compressed(destination/'writer.npz', scalars=np.asarray(writer_scalars),
        values=np.asarray(writer_vectors),offsets=offsets,scalar_names=WRITER_SCALARS,vector_names=WRITER_VECTORS,
        seeds=np.asarray(seed_rows),seed_fields=['entry',*SEED_FIELDS])
    if observed:
        np.savez_compressed(destination/'states.npz',values=np.asarray(states_out,dtype=np.float64),
            fields=STATE_FIELDS,geometry_ids=np.asarray(geometry_ids,dtype=np.uint64),
            hit_ids=np.asarray(hit_ids,dtype=np.uint64),tracks=np.asarray(tracks_out),track_fields=TRACK_FIELDS)
    summary=dict(source_events=task['events'],writer_tracks=ordinal,seed_count=len(seed_rows),
        audit_states=len(states_out),audit_observed=observed,
        first_source_entry=task['skip'],last_source_entry=task['skip']+len(seed_events)-1,
        seed_writer_event_entries=len(seed_events), reduction_version='uproot-v5',
        truth='Native writer SiHit-linked HepMC production p; explicit production only for operation diagnostics')
    if task.get('macro') and observed:
        (directory/'reference_3station_forward.root').symlink_to('track_states.root')
        invocation=task['macro']+'("reference_3station_forward.root","native_counts","afterGeo_latest","t_eQOP",false)'
        subprocess.run(['root','-l','-b','-q',invocation],check=True,cwd=directory)
        shutil = __import__('shutil')
        shutil.copy2(directory/'native_counts.root',destination/'native_counts.root')
    (destination/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
