"""Reduce actual Calypso audit states without storing the full JSON matrices."""
import json
import gzip
import shutil
from pathlib import Path

import numpy as np


STAGES = {'seed': 1, 'selected': 2, 'predicted': 3, 'filtered': 4, 'smoothed': 5,
          'converted': 6, 'refit_input': 7, 'refit_predicted': 8,
          'refit_filtered': 9, 'refit_smoothed': 10, 'final': 11}


def events(path):
    current = []
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else Path(path).open()) as source:
        for line in source:
            row = json.loads(line)
            if row['stage'] == 'event' and current:
                yield current
                current = []
            current.append(row)
    if current:
        yield current


def extract(path, output, task):
    rows, inputs, summaries, features = [], [], [], []
    def append(entry, record, code, barcode, kind, values, covariance=None, surface=None):
        values = np.asarray(values, dtype=float)
        if not np.all(np.isfinite(values)):
            return
        geometry = record.get('geometry_id', 0)
        if kind == 1:
            if values[4] == 0:
                return
            charge = np.sign(values[4])
        else:
            charge = record.get('charge', np.sign(record.get('bound', [0, 0, 0, 0, 1])[4]))
        layer = record.get('state', -1)
        identity = [entry, record['seed'], record['candidate'], barcode, code, layer, kind, geometry, charge]
        inputs.append(' '.join(map(str, identity))+ ' '+ ' '.join(format(x, '.17g') for x in values))
        rows.append(dict(entry=entry, seed=record['seed'], candidate=record['candidate'], barcode=barcode,
                         stage=code, state=layer, geometry=geometry, values=values, kind=kind, charge=charge,
                         covariance=np.asarray(covariance, dtype=float) if covariance is not None else np.full((6, 6), np.nan),
                         surface=np.asarray(surface, dtype=float) if surface is not None else np.full((4, 4), np.nan),
                         measurement=record.get('measurement_flag', False),
                         material_X0=record.get('material_full_X0', 0.),
                         material_path=record.get('material_path_correction', 1.),
                         path_mm=record.get('path_length_mm', np.nan),
                         update_error=np.asarray(record.get('_update_error', [np.nan]*6)),
                         smoothing_error=np.asarray(record.get('_smoothing_error', [np.nan]*6))))
    for ordinal, records in enumerate(events(path)):
        entry = task['skip'] + ordinal
        selected = {(r['seed'], r['candidate']): r for r in records if r['stage'] == 'selected'}
        majority = {key: r['matches'][0][0] if r.get('matches') else 0 for key, r in selected.items()}
        seed_barcode = {r['seed']: r['matches'][0][0] if r.get('matches') else 0 for r in records if r['stage'] == 'seed'}
        groups = {}
        for r in records:
            if r['stage'] in ['selected_state','refit_state']:
                groups.setdefault((r['stage'],r['seed'],r['candidate']),[]).append(r)
        for key,group in groups.items():
            ordered = list(reversed(group))
            for i,r in enumerate(ordered):
                if all(k in r for k in ['predicted','filtered','predicted_cov','projector','measurement','measurement_cov']):
                    p,f,c,h = (np.asarray(r[k],dtype=float) for k in ['predicted','filtered','predicted_cov','projector'])
                    innovation = h@c@h.T+np.asarray(r['measurement_cov'])
                    gain = np.linalg.solve(innovation,h@c).T
                    delta = f-p
                    delta[2] = (delta[2]+np.pi)%(2*np.pi)-np.pi
                    r['_update_error'] = delta-gain@(np.asarray(r['measurement'])-h@p)
                if i+1<len(ordered) and all(k in r for k in ['filtered','smoothed','filtered_cov']):
                    nxt = ordered[i+1]
                    if all(k in nxt for k in ['predicted_cov','jacobian','smoothed','predicted']):
                        cf,cp,j = (np.asarray(x,dtype=float) for x in [r['filtered_cov'],nxt['predicted_cov'],nxt['jacobian']])
                        try:
                            gain = np.linalg.solve(cp,j@cf).T
                            delta = np.asarray(nxt['smoothed'],dtype=float)-np.asarray(nxt['predicted'],dtype=float)
                            delta[2] = (delta[2]+np.pi)%(2*np.pi)-np.pi
                            difference = np.asarray(r['smoothed'],dtype=float)-np.asarray(r['filtered'],dtype=float)
                            difference[2] = (difference[2]+np.pi)%(2*np.pi)-np.pi
                            r['_smoothing_error'] = difference-gain@delta
                        except np.linalg.LinAlgError:
                            pass
        finals = {(r['seed'],r['candidate']):r for r in records if r['stage']=='final'}
        for key,r in selected.items():
            states = groups.get(('selected_state',*key),[])
            unique = {x['geometry_id']:x for x in states if 'material_full_X0' in x}
            paths = [x['path_length_mm'] for x in states if x.get('path_length_mm') is not None]
            seed_circle = next((x for x in records if x['stage']=='seed_circle' and x['seed']==key[0]),{})
            own_hits = set(r['hits'])
            other_hits = set().union(*(set(x['hits']) for k,x in selected.items() if k!=key))
            features.append(dict(entry=entry,seed=key[0],candidate=key[1],barcode=majority[key],
                purity=r['matches'][0][1]/max(1,len(r['hits'])) if r.get('matches') else 0.,
                measurements=r['measurements'],holes=r['holes'],outliers=r['outliers'],
                archived_path_span_mm=max(paths)-min(paths) if paths else np.nan,
                archived_full_slab_X0=sum(x['material_full_X0']*x.get('material_path_correction',1.) for x in unique.values()),
                archived_material_surfaces=len(unique),seed_stations=seed_circle.get('stations',-1),
                refit=finals.get(key,{}).get('refit',-1),shared_hit_fraction=len(own_hits&other_hits)/max(1,len(own_hits))))
        stage_counts = {}
        for r in records:
            name = r['stage']
            barcode = majority.get((r['seed'], r['candidate']), seed_barcode.get(r['seed'], 0))
            if name in ['seed', 'selected', 'converted', 'refit_input', 'final'] and 'momentum_A_GeV' in r:
                append(entry, r, STAGES[name], barcode, 0, [*r['position_A_mm'], *r['momentum_A_GeV']], r.get('cov'))
                stage_counts[name] = stage_counts.get(name, 0)+1
            if name not in ['selected_state', 'refit_state']:
                continue
            if not r.get('measurement_flag'):
                continue
            for component in ['predicted', 'filtered', 'smoothed']:
                if component in r:
                    code = STAGES[('refit_' if name == 'refit_state' else '')+component]
                    append(entry, r, code, barcode, 1, r[component], r.get(component+'_cov'), r.get('surface'))
        summaries.append(dict(source_entry=entry, selected=len(selected), stages=stage_counts))
    assert len(summaries) == task['events']
    output = Path(output)
    with Path('stage_states.npz').open('wb') as target:
        np.savez_compressed(target, **{key: np.asarray([r[key] for r in rows]) for key in rows[0]})
    np.savez_compressed('stage_features.npz',**{key:np.asarray([r[key] for r in features]) for key in features[0]})
    shutil.copy2('stage_states.npz',output/'states.npz')
    shutil.copy2('stage_features.npz',output/'candidate_features.npz')
    (output/'stage_inputs.txt').write_text('\n'.join(inputs)+'\n')
    (output/'stage_schema.json').write_text(json.dumps(dict(stages=STAGES, input_rows=len(rows),
        source_events=len(summaries), event_counts=summaries,
        reference='Native audit; no reset, covariance correction, truth seed or smoothing change unless explicitly selected',
        projection='InputStates passed back to the actual Calypso extrapolator',
        raw_JSON_retained=False), indent=2)+'\n')
    return len(rows)
