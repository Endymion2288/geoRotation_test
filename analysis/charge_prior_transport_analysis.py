#!/usr/bin/env python3
"""Compare a finite common physical prior through actual Calypso transport."""
import argparse
import json
from pathlib import Path

import numpy as np

from slope_tables import histogram, write_json


TARGETS = ['Common plane z_F=0 mm', 'VetoNu z_F=-3112 mm', 'First station-1 SiHit plane']


def quantiles(values):
    return dict(count=len(values), probabilities=[0., .5, .9, .99, 1.],
        values=np.quantile(values, [0., .5, .9, .99, 1.]).tolist() if len(values) else None)


def compare(a, b, cloud):
    for name in ['source_entry', 'seed', 'candidate', 'barcode', 'stage', 'layer', 'geometry_id', 'charge']:
        np.testing.assert_array_equal(a[name], b[name], err_msg=name)
    for data in [a, b]:
        np.testing.assert_array_equal(data['source_entry'], cloud['source_entry'])
        np.testing.assert_array_equal(data['candidate'], cloud['candidate'])
        np.testing.assert_allclose(data['input_position_F_mm'], cloud['position_F_mm'], rtol=0, atol=1e-10)
        np.testing.assert_allclose(data['input_momentum_F_GeV'], cloud['momentum_F_GeV'], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(data['input_qop_GeVinv'], cloud['qop_GeVinv'], rtol=1e-12, atol=1e-16)
        assert np.all(data['input_has_covariance'] == 0)
    rows = []
    for target, name in enumerate(TARGETS):
        sa, sb = a['projection_valid'][:, target], b['projection_valid'][:, target]
        common = (sa == 1) & (sb == 1)
        rows.append(dict(target=name, input_points=len(sa),
            status_codes=[-1, 0, 1],
            status_migration=np.bincount(3*(sa+1)+sb+1, minlength=9).reshape(3, 3).tolist(),
            before_success=int((sa == 1).sum()), after_success=int((sb == 1).sum()),
            common_success=int(common.sum())))
        signs = {}
        for label, data, good in [('beforeGeo', a, sa == 1), ('afterGeo', b, sb == 1)]:
            qop = data['projection_qop_GeVinv'][good, target]
            charge = data['projection_charge'][good, target]
            assert np.isfinite(qop).all() and np.all(qop != 0)
            np.testing.assert_array_equal(np.sign(qop), charge)
            signs[label] = dict(valid=len(qop),
                input_sign_changed=int(np.count_nonzero(charge != np.sign(data['input_qop_GeVinv'][good]))))
            p = np.linalg.norm(data['projection_F'][good, target, 3:], axis=1)
            np.testing.assert_allclose(np.abs(qop)*p, 1., rtol=1e-12, atol=1e-12)
        qa, qb = [x['projection_qop_GeVinv'][common, target] for x in [a, b]]
        signs['paired_endpoint_sign_differences'] = int(np.count_nonzero(np.sign(qa) != np.sign(qb)))
        va, vb = [x['projection_F'][common, target] for x in [a, b]]
        pa, pb = [np.linalg.norm(x[:, 3:], axis=1) for x in [va, vb]]
        ua, ub = [x[:, 3:]/p[:, None] for x, p in [(va, pa), (vb, pb)]]
        difference = dict(position_mm=np.linalg.norm(vb[:, :3]-va[:, :3], axis=1),
            direction_rad=np.arctan2(np.linalg.norm(np.cross(ua, ub), axis=1), np.einsum('ij,ij->i', ua, ub)),
            relative_p=pb/pa-1., relative_qop=qb/qa-1.)
        rows[-1].update(sign=signs, metrics={key: quantiles(np.abs(v)) for key, v in difference.items()},
            histograms={key: histogram(v, edges) for key, v, edges in [
                ('position_mm', difference['position_mm'], np.geomspace(1e-14, 1e-3, 111)),
                ('direction_rad', difference['direction_rad'], np.geomspace(1e-17, 1e-6, 111)),
                ('relative_p', difference['relative_p'], np.linspace(-.01, .01, 201)),
                ('relative_qop', difference['relative_qop'], np.linspace(-.01, .01, 201))]})
        extreme = np.argsort(np.abs(difference['relative_p']))[-5:][::-1]
        selected = np.flatnonzero(common)[extreme]
        rows[-1]['largest_relative_p_differences'] = [dict(candidate=int(cloud['candidate'][i]),
            source_entry=int(cloud['source_entry'][i]), input_qop_GeVinv=float(cloud['qop_GeVinv'][i]),
            input_p_GeV=float(np.linalg.norm(cloud['momentum_F_GeV'][i])),
            relative_p=float(difference['relative_p'][j]),
            position_difference_mm=float(difference['position_mm'][j]),
            direction_difference_rad=float(difference['direction_rad'][j])) for j, i in zip(extreme, selected)]
    return rows


def reduce(source, output):
    manifest = json.loads((source/'manifest.json').read_text())
    previous = Path(manifest['previous'])
    groups, validation, total_sign_changes = {}, {}, 0
    for dsid in manifest['samples']:
        with np.load(previous/dsid/'physical_cloud.npz') as file:
            cloud = dict(file)
        n = len(cloud['candidate'])
        assert len(np.unique(cloud['candidate'])) == n
        for control in manifest['controls']:
            tables = []
            for variant in ['beforeGeo', 'afterGeo']:
                path = source/dsid/control/variant
                assert (path/'exit_code.txt').read_text().strip() == '0'
                run = json.loads((path/'run.json').read_text())
                assert run['success'] and run['args']['material'] == manifest['materials'][variant]
                assert run['interaction_energy_loss'] == (control == 'material_on')
                assert run['interaction_scattering'] == (control == 'material_on')
                actual = [Path(x) for x in run['loaded_acts_libraries']]
                expected = Path(manifest['endpoint_observer']).parent/'lib/libCalypsoSlopeDiagnostics.so'
                assert expected.resolve() in [p.resolve() for p in actual]
                for name in ['libFaserActsGeometry.so', 'libFaserActsGeometryLib.so']:
                    expected_geometry = Path(manifest['libraries'][variant])/name
                    assert expected_geometry.resolve() in [p.resolve() for p in actual]
                with np.load(path/'transport.npz') as file:
                    data = {k: file[k][np.argsort(file['candidate'])] for k in file.files}
                assert len(data['candidate']) == n
                fields = []
                if control == 'material_on':
                    with np.load(previous/dsid/variant/'transport.npz') as old:
                        order = np.argsort(old['candidate'])
                        for name in old.files:
                            np.testing.assert_array_equal(old[name][order], data[name], err_msg=dsid+'/'+variant+'/'+name)
                            fields.append(name)
                validation['|'.join([dsid, control, variant])] = dict(actual_library_paths=[str(p) for p in actual],
                    original_fields_exact=fields, source=str(path), input_points=n)
                tables.append(data)
            rows = compare(*tables, cloud)
            for row in rows:
                total_sign_changes += row['sign']['paired_endpoint_sign_differences']
                total_sign_changes += sum(row['sign'][v]['input_sign_changed'] for v in ['beforeGeo', 'afterGeo'])
            groups[dsid+'|'+control] = rows
    output.mkdir(exist_ok=False)
    write_json(output/'summary.json', dict(groups=groups, validation=validation,
        total_observed_sign_changes=total_sign_changes, source_manifest=str(source/'manifest.json'),
        input_points_per_species=manifest['samples_per_seed']*5,
        definition='Same full finite physical kinematic prior draws on z_F=0 mm, exact nonlinear rotation, native Calypso extrapolators and current/strict-rotated material pair',
        scope='Transport only. Five representative seed priors per species, not a reconstruction efficiency or misidentification-rate measurement. The time marginal is integrated out for this static transport. No stochastic G4 ensemble or CKF covariance equivalence is asserted.'))
    print('Actual endpoint sign changes:', total_sign_changes)


def plot(source, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data = json.loads(source.read_text())
    controls = data.get('plot_controls', [
        dict(key='material_on', label='material on', color='#c43c39'),
        dict(key='material_off', label='material off', color='#205fb2')])
    output.mkdir(exist_ok=False)
    index, counts = [], {}
    fields = [('position_mm', 'Position difference [mm]', True),
        ('direction_rad', 'Direction difference [rad]', True),
        ('relative_p', r'$(p_{after}-p_{before})/p_{before}$', False)]
    with PdfPages(output/'physical_prior_transport_report.pdf') as pdf:
        for dsid in ['100044', '100116', '100117']:
            fig, axes = plt.subplots(3, 3, figsize=(13, 9))
            for target in range(3):
                for ax, (field, xlabel, log) in zip(axes[target], fields):
                    flows = []
                    for definition in controls:
                        control, color = definition['key'], definition['color']
                        h = data['groups'][dsid+'|'+control][target]['histograms'][field]
                        assert sum(h['counts'])+h['underflow']+h['overflow'] == h['finite']
                        ax.stairs(h['counts'], h['edges'], color=color, label=definition['label'])
                        flows.append(h['underflow'])
                        counts[f'{dsid}|{control}|{target}|{field}'] = h
                    if log:
                        ax.set_xscale('log')
                        ax.set_ylim(bottom=0)
                    else:
                        ax.set_yscale('log')
                    flow_label = 'Underflow on/off' if 'plot_controls' not in data else 'Underflow legacy/row-major'
                    below_legend = field == 'position_mm' and target == 0
                    ax.text(.98, .03 if below_legend else .98, f'{flow_label}: {flows[0]}/{flows[1]}',
                        transform=ax.transAxes, ha='right', va='bottom' if below_legend else 'top', fontsize=7)
                    if field == 'position_mm' and target == 0:
                        ax.set_ylim(0, 1)
                        ax.set_yticks([0, 1])
                        ax.text(.5, .5, 'All position differences are exactly zero',
                            transform=ax.transAxes, ha='center', fontsize=8)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel('Entries / bin')
                    ax.set_title(TARGETS[target], fontsize=9)
                    ax.grid(alpha=.2)
            axes[0, 0].legend(fontsize=7, loc='upper left')
            title = data.get('plot_title', 'native Calypso before/after rotation')
            fig.suptitle(f'{dsid}: identical finite physical prior | {title} | 2,560 points', fontsize=12)
            fig.text(.5, .025, '1 RDO file, source entries [0,20) | 5 representative seed priors, 512 draws each | time marginalized',
                ha='center', fontsize=8)
            fig.text(.5, .01, 'Transport only | actual endpoint charge and q/p recorded', ha='center', fontsize=8)
            fig.tight_layout(rect=[0, .055, 1, .95])
            suffix = 'material_decoder_transport' if 'plot_controls' in data else 'physical_prior_transport'
            name = dsid+'_'+suffix
            fig.savefig(output/(name+'.png'), dpi=150)
            fig.savefig(output/(name+'.pdf'))
            pdf.savefig(fig)
            plt.close(fig)
            index.append(dict(name=name, png=name+'.png', pdf=name+'.pdf',
                reference='S06 Q5 finite common physical prior; '+title))
    write_json(output/'figure_index.json', index)
    write_json(output/'bin_counts.json', counts)
    write_json(output/'plot_source.json', dict(summary=str(source), pages=len(index)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--phase', choices=['reduce', 'plot'], required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (reduce if args.phase == 'reduce' else plot)(args.input, args.output)


if __name__ == '__main__':
    main()
