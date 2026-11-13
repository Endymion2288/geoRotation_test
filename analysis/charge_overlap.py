#!/usr/bin/env python3
"""Truth-only overlap weighting, learned on discovery and tested on held-out files."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from charge_tables import weighted_rate
from slope_tables import write_json


DSIDS = ['100116', '100117']
ENTRANCE_FEATURES = ['log production p', 'st1 truth x', 'st1 truth y', 'sensor true sx', 'sensor true sy']
PRODUCTION_FEATURES = ['log production p', 'production x', 'production y', 'production sx', 'production sy']


def eligible(data):
    return data.get('population', data['fiducial']) & np.isfinite(data['truth_features']).all(axis=1)


def basis(x, mean, scale, knots):
    z = (x-mean)/scale
    columns = [np.ones(len(z))]
    columns.extend(z[:, j] for j in range(5))
    columns.extend(z[:, j]**2 for j in range(5))
    columns.extend(z[:, i]*z[:, j] for i in range(5) for j in range(i+1, 5))
    columns.extend(np.maximum(z[:, j]-k, 0) for j in range(5) for k in knots[j])
    return np.column_stack(columns)


def fit_model(data, features=ENTRANCE_FEATURES):
    masks = {k: eligible(data[k]) &
             (data[k]['file_index'] < 3) for k in DSIDS}
    x = np.concatenate([data[k]['truth_features'][masks[k]] for k in DSIDS])
    mean, scale = x.mean(axis=0), x.std(axis=0)
    assert (scale > 0).all()
    knots = np.quantile((x-mean)/scale, [.2, .5, .8], axis=0).T
    design = basis(x, mean, scale, knots)
    n = [int(masks[k].sum()) for k in DSIDS]
    y = np.r_[np.ones(n[0]), np.zeros(n[1])]
    class_weight = np.r_[np.full(n[0], len(y)/(2*n[0])), np.full(n[1], len(y)/(2*n[1]))]
    penalty = 1e-6

    def objective(beta):
        logit = design@beta
        likelihood = np.mean(class_weight*(np.logaddexp(0, logit)-y*logit))
        gradient = design.T@(class_weight*(expit(logit)-y))/len(y)
        regularizer = beta.copy()
        regularizer[0] = 0
        return likelihood+.5*penalty*(regularizer@regularizer), gradient+penalty*regularizer

    fit = minimize(objective, np.zeros(design.shape[1]), jac=True, method='L-BFGS-B',
                   options=dict(maxiter=1000, ftol=1e-12, gtol=1e-8))
    report = dict(coefficients=fit.x.tolist(), mean=mean.tolist(), scale=scale.tolist(),
        knots=knots.tolist(), basis='Intercept, five standardized inputs, squares, all pairwise products, hinges at discovery 20/50/80 percentiles',
        penalty=penalty, optimizer_success=bool(fit.success), optimizer_message=str(fit.message),
        optimizer_iterations=int(fit.nit), loss=float(fit.fun),
        max_gradient=float(np.max(np.abs(fit.jac))), training_truth_per_charge=dict(zip(DSIDS, n)),
        features=features,
        learning_inputs='Truth only; no reconstructed p, charge estimate, quality or correctness',
        class_prior='Equal training charge weights',
        target='Overlap density proportional to f_plus(x)*f_minus(x)/(f_plus(x)+f_minus(x))',
        probability_support=[.01, .99], held_out_files=[3, 4])
    assert fit.success, report
    return report


def weighted_ks(a, wa, b, wb):
    order_a, order_b = np.argsort(a), np.argsort(b)
    a, b, wa, wb = a[order_a], b[order_b], wa[order_a], wb[order_b]
    grid = np.unique(np.r_[a, b])
    ca, cb = np.r_[0., np.cumsum(wa)/wa.sum()], np.r_[0., np.cumsum(wb)/wb.sum()]
    return float(np.max(np.abs(ca[np.searchsorted(a, grid, side='right')]-
                              cb[np.searchsorted(b, grid, side='right')])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--feature-file', type=Path, help='Production-only features with exact source identity')
    parser.add_argument('--population', choices=['fiducial', 'all-primary'], default='fiducial')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    data = {}
    for k in DSIDS:
        with np.load(args.input/f'{k}_afterGeo_outcomes.npz') as archive:
            data[k] = {name: archive[name] for name in archive.files}
        data[k]['population'] = data[k]['fiducial'] if args.population == 'fiducial' else np.ones(len(data[k]['fiducial']), dtype=bool)
        if args.feature_file:
            with np.load(args.feature_file) as features:
                np.testing.assert_array_equal(data[k]['sample_entry'], features[k+'_sample_entry'])
                data[k]['truth_features'] = features[k]
    model = fit_model(data, PRODUCTION_FEATURES if args.feature_file else ENTRANCE_FEATURES)
    model['population'] = args.population
    model['feature_file'] = str(args.feature_file) if args.feature_file else None
    write_json(args.output/'model.json', model)
    weights, propensity = {}, {}
    mean, scale, knots = (np.array(model[k]) for k in ['mean', 'scale', 'knots'])
    for k in DSIDS:
        valid = eligible(data[k])
        prob = np.full(len(valid), np.nan)
        prob[valid] = expit(basis(data[k]['truth_features'][valid], mean, scale, knots)@np.array(model['coefficients']))
        propensity[k] = prob
        supported = valid & (prob >= .01) & (prob <= .99)
        weights[k] = np.where(supported, 1-prob if k == '100116' else prob, 0.)
    balance, samples, rates = {}, {}, {}
    for split in ['discovery', 'held_out', 'all']:
        masks = {k: eligible(data[k]) &
            (np.ones(len(weights[k]), dtype=bool) if split == 'all' else
             data[k]['file_index'] < 3 if split == 'discovery' else data[k]['file_index'] >= 3) for k in DSIDS}
        smd, ks, original = [], [], []
        for j in range(5):
            x = [data[k]['truth_features'][masks[k], j] for k in DSIDS]
            w = [weights[k][masks[k]] for k in DSIDS]
            scale_j = np.sqrt(.5*(np.var(x[0])+np.var(x[1])))
            original.append(float((x[0].mean()-x[1].mean())/scale_j))
            smd.append(float((np.average(x[0], weights=w[0])-np.average(x[1], weights=w[1]))/scale_j))
            ks.append(weighted_ks(x[0], w[0], x[1], w[1]))
        balance[split] = dict(original_SMD=original, overlap_SMD=smd, marginal_weighted_KS_distance=ks,
            passes_frozen_balance_threshold=max(np.abs(smd)) < .05 and max(ks) < .05,
            criterion='All five |SMD|<0.05 and weighted marginal KS distances<0.05')
        samples[split] = {}
        rates[split] = {}
        for k in DSIDS:
            w = weights[k]*masks[k]
            positive = w > 0
            samples[split][k] = dict(eligible=int(masks[k].sum()), supported=int(positive.sum()),
                sumw=float(w.sum()), sumw2=float(w@w), effective_n=float(w.sum()**2/(w@w)),
                weight_quantiles=np.quantile(w[positive], [.5, .9, .95, .99, 1]).tolist())
        for variant in ['beforeGeo', 'afterGeo']:
            rates[split][variant] = {}
            for k in DSIDS:
                with np.load(args.input/f'{k}_{variant}_outcomes.npz') as archive:
                    selected, wrong, status = archive['selected'], archive['wrong'], archive['status']
                w = weights[k]*masks[k]
                rates[split][variant][k] = dict(conditional=weighted_rate(wrong, selected, w),
                    truth_correct_efficiency=weighted_rate(status == 1, data[k]['fiducial'], w),
                    truth_wrong_efficiency=weighted_rate(wrong, data[k]['fiducial'], w),
                    fiducial_acceptance=weighted_rate(data[k]['fiducial'], masks[k], w),
                    full_population_correct_efficiency=weighted_rate(status == 1, masks[k], w))
    differences = {}
    for split in rates:
        differences[split] = {}
        for variant in rates[split]:
            a, b = [rates[split][variant][k]['conditional'] for k in ['100117', '100116']]
            delta = a['value']-b['value']
            error = float(np.hypot(a['standard_error'], b['standard_error']))
            differences[split][variant] = dict(mu_minus_minus_mu_plus=delta,
                standard_error=error, interval95=[delta-1.96*error, delta+1.96*error])
    write_json(args.output/'summary.json', dict(model=model, balance=balance, samples=samples,
        rates=rates, differences=differences,
        interval_scope='Source-event outcome uncertainty conditional on the discovery-frozen propensity model; training-model uncertainty is not included'))
    np.savez_compressed(args.output/'truth_weights.npz', **{k: weights[k] for k in DSIDS})
    print(json.dumps(dict(balance=balance, differences=differences), indent=2), flush=True)


if __name__ == '__main__':
    main()
