"""Reduce Calypso's existing ROOT writers; no fitted parameters are recalculated."""
import json
from pathlib import Path
import shutil

import awkward as ak
import numpy as np
import uproot


LABELS = {1: 'Seed population', 3: 'CKF predicted', 4: 'CKF filtered',
          5: 'CKF smoothed', 8: 'Refit predicted', 9: 'Refit filtered',
          10: 'Refit smoothed', 11: 'Final'}


def tree_arrays(path, branches):
    with uproot.open(path) as file:
        key = next(k for k, v in file.classnames().items() if v == 'TTree')
        return file[key].arrays(branches, library='ak')


def extract(output, task, arrays):
    events, tracks = arrays['events'], arrays['tracks']
    # Production IDs repeat across units, but the frozen diagnostic shards never
    # cross a repeat. CKF row order is additionally checked against output counts.
    event_ids = events['event'].tolist()
    assert len(set(event_ids)) == len(event_ids), 'Split diagnostics at production-ID repeats'
    entries = dict(zip(event_ids, events['source_entry'].tolist()))
    rows, inputs, angle_rows = [], [], []
    checks = dict(source_events=len(event_ids), unmatched_refit_rows=0, ambiguous_refit_rows=0)
    rotated = task['variant'] != 'beforeGeo'

    def append(entry, ordinal, barcode, code, values, charge, station=-1, layer=-1):
        values = np.asarray(values, float)
        if not np.all(np.isfinite(values)) or np.linalg.norm(values[3:]) == 0:
            return
        rows.append(dict(entry=entry, track_index=ordinal, barcode=barcode, stage=code,
                         values=values, charge=charge, station=station, layer=layer))
        identity = [entry, ordinal, 0, barcode, code, layer, 0, 0, charge]
        inputs.append(' '.join(map(str, identity)) + ' ' + ' '.join(format(x, '.17g') for x in values))

    def signature(data, row):
        keys = ['station', 'layer', 'eta_module', 'phi_module', 'side', 'l_x_hit']
        cols = [ak.to_numpy(data[k][row]) for k in keys]
        return tuple(sorted(tuple(float(x) if j < 5 else round(float(x), 9)
                                  for j, x in enumerate(v)) for v in zip(*cols)))

    parts = {}
    base = ['event_nr', 'track_nr', 'majorityParticleId', 'station', 'layer',
            'eta_module', 'phi_module', 'side', 'l_x_hit']
    for suffix in ['prt', 'flt', 'smt']:
        base += [f'{name}_{suffix}' for name in ['g_x', 'g_y', 'g_z', 'px', 'py', 'pz', 'eQOP']]
        base += [f'{name}_{suffix}' for name in ['ePHI','eTHETA','err_ePHI','err_eTHETA','err_eQOP']]
    for label, filename in [('ckf', 'slope_3station_forward_track_states_ckf.root'),
                            ('refit', 'slope_3station_forward_track_states_ckf1.root')]:
        parts[label] = tree_arrays(filename, base)
    ckf = parts['ckf']
    expected = np.repeat(events['event'], events['n_tracks'])
    assert np.array_equal(ak.to_numpy(ckf['event_nr']), expected)
    assert len(ckf) == len(tracks['source_entry'])
    assert np.array_equal(np.array([entries[int(x)] for x in ak.to_list(ckf['event_nr'])]), tracks['source_entry'])
    lookup = {}
    for row in range(len(ckf)):
        key = (int(ckf['event_nr'][row]), signature(ckf, row))
        lookup.setdefault(key, []).append(row)
    for label, data in parts.items():
        checks[label + '_tracks'] = len(data)
        for row in range(len(data)):
            event_id = int(data['event_nr'][row])
            if label == 'ckf':
                track_row = row
            else:
                matches = lookup.get((event_id, signature(data, row)), [])
                if len(matches) != 1:
                    checks['unmatched_refit_rows' if not matches else 'ambiguous_refit_rows'] += 1
                    continue
                track_row = matches[0]
            # Final collection preserves selected CKF order, including refit
            # fallback. The final majority remains fixed across P/F/S comparisons.
            entry = int(tracks['source_entry'][track_row])
            ordinal = int(tracks['track_index'][track_row]) if 'track_index' in tracks else int(ckf['track_nr'][track_row])
            barcode = int(tracks['barcode'][track_row])
            # These are original float writer parameters, before any additional
            # transport or angle normalization in this analysis.
            angle_columns={field+'_'+suffix:ak.to_numpy(data[field+'_'+suffix][row])
                for suffix in ['prt','flt','smt']
                for field in ['ePHI','eTHETA','eQOP','err_ePHI','err_eTHETA','err_eQOP']}
            modules={field:ak.to_numpy(data[field][row]) for field in ['station','layer','side']}
            for hit in range(len(modules['station'])):
                record=dict(entry=entry,track_index=ordinal,barcode=barcode,
                            writer=0 if label=='ckf' else 1,
                            **{field:int(values[hit]) for field,values in modules.items()})
                record.update({field:float(values[hit]) for field,values in angle_columns.items()})
                angle_rows.append(record)
            for j, suffix in enumerate(['prt', 'flt', 'smt']):
                values = np.column_stack([ak.to_numpy(data[f'{k}_{suffix}'][row])
                                          for k in ['g_x', 'g_y', 'g_z', 'px', 'py', 'pz']])
                finite = np.all(np.isfinite(values), axis=1)
                if not finite.any():
                    continue
                # The same most-upstream measurement convention as the audit
                # analysis; the actual Calypso tool transports it to fixed planes.
                physical_z = values[:, 0] if rotated else values[:, 2]
                index = np.flatnonzero(finite)[np.argmin(physical_z[finite])]
                append(entry, ordinal, barcode, (3 if label == 'ckf' else 8) + j,
                       values[index], np.sign(float(data[f'eQOP_{suffix}'][row][index])),
                       int(data['station'][row][index]), int(data['layer'][row][index]))
    for row in range(len(tracks['source_entry'])):
        ordinal = int(tracks['track_index'][row]) if 'track_index' in tracks else int(ckf['track_nr'][row])
        append(int(tracks['source_entry'][row]), ordinal, int(tracks['barcode'][row]),
               11, tracks['upstream_raw'][row], float(tracks['charge'][row]))
    seed = tree_arrays('slope_3station_forward_seed_summary_circleFitTrackSeedTool.root',
        ['event_number', 'majorityParticleId', 'track_vx', 'track_vy', 'track_vz',
         'track_px', 'track_py', 'track_pz', 'track_charge'])
    for row in range(len(seed)):
        entry = entries[int(seed['event_number'][row])]
        for i in range(len(seed['track_px'][row])):
            append(entry, -(i+1), int(seed['majorityParticleId'][row][i]), 1,
                   [float(seed[k][row][i]) for k in ['track_vx', 'track_vy', 'track_vz', 'track_px', 'track_py', 'track_pz']],
                   float(seed['track_charge'][row][i]))
    np.savez_compressed('native_states.npz', **{k: np.asarray([r[k] for r in rows]) for k in rows[0]})
    shutil.copy2('native_states.npz', output/'native_states.npz')
    np.savez_compressed('native_angles.npz', **{k: np.asarray([r[k] for r in angle_rows]) for k in angle_rows[0]})
    shutil.copy2('native_angles.npz', output/'native_angles.npz')
    (output/'native_inputs.txt').write_text('\n'.join(inputs) + '\n')
    checks.update(stage_labels=LABELS, rows=len(rows), seed_association='Population only; not a same-seed claim',
                  refit_matching='Same source event and complete module/local measurement signature',
                  units='Raw native frame; mm, GeV; independent Calypso transport required')
    (output/'native_state_validation.json').write_text(json.dumps(checks, indent=2) + '\n')
    return len(rows)
