#!/usr/bin/env python3
"""Replay direction-chart quality metrics from saved actual Calypso outputs."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from slope_centroids import collect
from slope_periodic_refit import CONTROLS
from slope_tables import join_truth, write_json
from validate_slope_refit_artifacts import histograms, pages, page_text


def check_histogram(values, saved):
    finite = values[np.isfinite(values)]
    edges = np.asarray(saved['edges'])
    bins = np.searchsorted(edges, finite, side='right') - 1
    inside = (bins >= 0) & (bins < len(edges)-1)
    counts = np.bincount(bins[inside], minlength=len(edges)-1)
    for name in ['counts', 'sumw', 'sumw2']:
        np.testing.assert_array_equal(counts, saved[name])
    assert saved['raw_entries'] == len(values)
    assert saved['finite'] == len(finite)
    assert saved['invalid'] == len(values)-len(finite)
    assert saved['underflow'] == int((bins < 0).sum())
    assert saved['overflow'] == int((bins >= len(edges)-1).sum())


def validate(summary):
    data = json.loads(summary.read_text())
    folder = summary.parent
    report = folder/'chart_diagnostics_report.pdf'
    assert pages(report) == 2*len(data['cases'])
    checks, page, histogram_count = {}, 0, 0
    for label, case in data['cases'].items():
        native, count = collect([Path(p) for p in case['native_paths']])
        control, other_count = collect([Path(p) for p in case['control_paths']])
        assert count == other_count == case['source_events']
        r, q, truth = native['tracks'], control['tracks'], native['truth']
        for key in ['source_entry', 'ordinal']:
            np.testing.assert_array_equal(r[key], q[key])
        same = (r['hit_signatures'] == q['hit_signatures']) & (r['barcode'] == q['barcode'])
        index, linked = join_truth(native)
        selected = linked & (np.abs(truth['pdg'][index]) == 13) & (r['purity'] >= .7) & ((r['station_bits'] & 14) == 14)
        if case['validation_source_entry_min'] is not None:
            selected &= r['source_entry'] >= case['validation_source_entry_min']
        assert int(selected.sum()) == case['selected_candidates']
        p = np.linalg.norm(r['upstream_F'][:, 3:], axis=1)
        pt = np.linalg.norm(truth['production_momentum_F_GeV'][index], axis=1)
        windows = dict(all_p=np.ones(len(p), bool), native_p=(p > 50) & (p < 250),
                       truth_p=(pt > 50) & (pt < 250))
        windows.update({k+'|same_hits_truth': v & same for k, v in list(windows.items())})
        mode = case['control']
        for plane, i in [('st1', 0), ('sensor', 2)]:
            with np.errstate(divide='ignore', invalid='ignore'):
                values = {'native': np.where(r['projection_valid'][:, i] == 1,
                    r['projection'][:, i, 4]/r['projection'][:, i, 5], np.nan),
                    mode: np.where(q['projection_valid'][:, i] == 1,
                    q['projection'][:, i, 4]/q['projection'][:, i, 5], np.nan)}
                reference = (truth['production_momentum_F_GeV'][index] if i == 0 and case['generated_st1_truth']
                    else truth['truth_model_momentum_F_GeV'][index, i] if i == 0
                    else truth['truth_sensor_direction_F'][index])
                values['truth'] = reference[:, 1]/reference[:, 2]
            finite = selected & np.logical_and.reduce([np.isfinite(v) for v in values.values()])
            item = case['planes'][plane]
            for name, window in windows.items():
                mask = finite & window
                saved = item['selections'][name]
                assert saved['common_finite'] == int(mask.sum())
                for width, row in saved['windows'].items():
                    actual = {k: int((mask & (np.abs(v) < int(width)*1.e-6)).sum()) for k, v in values.items()}
                    assert actual == row['counts']
                    for key, first, second in [('native_minus_truth', 'native', 'truth'),
                        ('control_minus_truth', mode, 'truth'), ('control_minus_native', mode, 'native')]:
                        interval = row[key]
                        numerator = actual[first]-actual[second]
                        assert interval['numerator'] == numerator and interval['denominator'] == int(mask.sum())
                        np.testing.assert_allclose(interval['estimate'], numerator/mask.sum(), rtol=0, atol=1.e-16)
                for key in ['native', mode]:
                    residual = values[key]-values['truth']
                    out = saved['residuals'][key]
                    check_histogram(residual[mask], out['histogram'])
                    sample = residual[mask]
                    np.testing.assert_allclose(out['summary']['rms'], np.linalg.norm(sample)/np.sqrt(len(sample)), rtol=1.e-14)
                    np.testing.assert_allclose(out['summary']['absolute_error_68'], np.quantile(np.abs(sample), .68), rtol=1.e-14)
                    source = r if key == 'native' else q
                    variance = source['projection'][:, i, 6]
                    valid_pull = mask & np.isfinite(variance) & (variance > 0)
                    check_histogram(residual[valid_pull]/np.sqrt(variance[valid_pull]), out['pull_histogram'])
                    assert out['invalid_pull_variance'] == int((mask & ~valid_pull).sum())
                    for width, tail in out['large_error'].items():
                        assert tail['numerator'] == int((mask & (np.abs(residual) > float(width)*.001)).sum())
                        assert tail['denominator'] == int(mask.sum())
            for title, key in [('Native', 'native'), ('Truth', 'truth'), (CONTROLS[mode][0], mode)]:
                check_histogram(values[key][finite & windows['native_p']], item['figure_histograms'][title])
            if case['after_paths']:
                after, after_count = collect([Path(p) for p in case['after_paths']])
                assert after_count == count
                ar, at = after['tracks'], after['truth']
                ai, found = join_truth(after)
                ap = np.linalg.norm(ar['upstream_F'][:, 3:], axis=1)
                chosen = found & (np.abs(at['pdg'][ai]) == 13) & (ar['purity'] >= .7) & ((ar['station_bits'] & 14) == 14)
                chosen &= (ap > 50) & (ap < 250) & (ar['projection_valid'][:, i] == 1)
                with np.errstate(divide='ignore', invalid='ignore'):
                    sy = ar['projection'][:, i, 4]/ar['projection'][:, i, 5]
                check_histogram(sy[chosen], item['figure_histograms']['afterGeo_own_cohort'])
                assert item['afterGeo']['selected'] == int(chosen.sum())
                assert item['afterGeo']['near_zero'] == int((chosen & (np.abs(sy) < .0001)).sum())
            stem = f'{label}_{plane}_{mode}_diagnostics'
            png, pdf = folder/(stem+'.png'), folder/(stem+'.pdf')
            with Image.open(png) as image:
                image.load()
            assert pages(pdf) == 1
            page += 1
            assert page_text(pdf) == page_text(report, page)
        number = sum(1 for _ in histograms(case))
        histogram_count += number
        checks[label] = dict(source_events=count, replayed_selections=12,
            replayed_window_counts=True, replayed_residuals_and_pulls=True,
            verified_histograms=number)
    assert len(list(folder.glob('*.png'))) == page
    return dict(summary=str(summary), cases=checks, pages=page,
        histograms=histogram_count, png_decoding=True, single_report_page_text_equal=True,
        definition='Raw saved Calypso arrays replay, independent searchsorted/bin accumulation, source identities and fixed native selection; no tracking or propagation recomputation')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.summary)
    write_json(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
