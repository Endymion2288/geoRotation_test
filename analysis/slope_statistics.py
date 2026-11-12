"""Event-level uncertainty and independent conditional responses for S05."""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks, peak_widths


def event_ratio_interval(entries, numerator, denominator, event_range, seed=50519, repeats=1000):
    """Exact empirical event bootstrap, compressed by repeated count vectors."""
    begin, end = event_range
    entries = np.asarray(entries)-begin
    keep = (entries >= 0) & (entries < end-begin)
    num = np.bincount(entries[keep], weights=np.asarray(numerator)[keep], minlength=end-begin)
    den = np.bincount(entries[keep], weights=np.asarray(denominator)[keep], minlength=end-begin)
    patterns, counts = np.unique(np.stack([num,den],axis=1), axis=0, return_counts=True)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(end-begin, counts/(end-begin), size=repeats) @ patterns
    good = draws[:,1] > 0
    estimate = float(num.sum()/den.sum()) if den.sum() else None
    intervals = np.quantile(draws[good,0]/draws[good,1],[.025,.975]).tolist() if good.any() else [None,None]
    return dict(numerator=float(num.sum()), denominator=float(den.sum()), estimate=estimate,
                interval_95=intervals, resampling_unit='source event', bootstrap_repeats=repeats, seed=seed)


def distribution_structure(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    unique, counts = np.unique(finite, return_counts=True)
    order = np.argsort(counts)[::-1][:15]
    edges = np.linspace(-.012,.012,2401)
    grid = .5*(edges[1:]+edges[:-1])
    fine = np.histogram(finite,edges)[0].astype(float)
    kde = {}
    for bandwidth in [5e-5,1e-4,2e-4,4e-4]:
        density = gaussian_filter1d(fine,bandwidth/(edges[1]-edges[0]),mode='constant')
        density /= max(1,len(finite))*(edges[1]-edges[0])
        kde[f'{bandwidth:g}'] = density.tolist()
    quantiles = np.linspace(0,1,min(20001,len(finite))) if len(finite) else np.array([])
    return dict(finite=len(finite), exact_unique=len(unique),
        exact_zero=int(np.count_nonzero(finite==0)),
        within_1e_minus_12_of_zero=int(np.count_nonzero(np.abs(finite)<1e-12)),
        most_repeated=[dict(value=float(unique[i]),count=int(counts[i])) for i in order],
        unique_spacing_quantiles=np.quantile(np.diff(unique),[0,.01,.5,.99,1]).tolist() if len(unique)>1 else [],
        ecdf=dict(x=np.quantile(finite,quantiles).tolist() if len(finite) else [],cdf=quantiles.tolist()),
        kde=dict(grid=grid.tolist(),curves=kde,method='Gaussian convolution of 1e-5-wide histogram; includes finite out-of-plot values'))


def discover_windows(values):
    """Freeze all prominent discovery modes; never select using validation."""
    edges = np.linspace(-.01,.01,1001)
    centers = .5*(edges[1:]+edges[:-1])
    counts = gaussian_filter1d(np.histogram(values[np.isfinite(values)],edges)[0].astype(float),5.)
    indices, info = find_peaks(counts,prominence=max(5.,.025*counts.max(initial=0)),distance=15)
    widths = peak_widths(counts,indices,rel_height=.5)[0]*(edges[1]-edges[0]) if len(indices) else []
    return [dict(center=float(centers[i]),half_width=float(np.clip(w/2,1e-4,5e-4)),
                 fwhm_discovery=float(w),prominence_discovery=float(prom))
            for i,w,prom in zip(indices,widths,info['prominences'])]


def peak_statistics(values, entries, selected, windows, event_range):
    values = np.asarray(values)
    selected = selected & np.isfinite(values)
    results = []
    for window in windows:
        center, width = window['center'],window['half_width']
        peak = selected & (np.abs(values-center)<width)
        side = selected & (np.abs(values-center)>=width) & (np.abs(values-center)<3*width)
        fraction = event_ratio_interval(entries,peak,selected,event_range)
        excess = event_ratio_interval(entries,peak.astype(float)-.5*side,selected,event_range)
        in_range = (entries>=event_range[0]) & (entries<event_range[1])
        sample = values[peak & in_range]
        # Cluster delta-method intervals retain the correlation of tracks in
        # one source event. The reported width is within the frozen window.
        moment_interval = dict(mean_95=[None,None],std_95=[None,None])
        if len(sample)>1:
            mean, variance = np.mean(sample),np.var(sample)
            ids = entries[peak & in_range]-event_range[0]
            nevents = event_range[1]-event_range[0]
            mean_influence = np.bincount(ids,weights=sample-mean,minlength=nevents)/len(sample)
            var_influence = np.bincount(ids,weights=(sample-mean)**2-variance,minlength=nevents)/len(sample)
            correction = nevents/max(1,nevents-1)
            mean_se = np.sqrt(correction*np.sum(mean_influence**2))
            std_se = np.sqrt(correction*np.sum(var_influence**2))/(2*np.sqrt(variance)) if variance else 0.
            moment_interval = dict(mean_95=[float(mean-1.96*mean_se),float(mean+1.96*mean_se)],
                std_95=[float(max(0,np.sqrt(variance)-1.96*std_se)),float(np.sqrt(variance)+1.96*std_se)],
                method='Source-event cluster delta method; conditional moments inside the frozen window')
        results.append(dict(window=window,fraction=fraction,sideband_subtracted_fraction=excess,
            peak_mean=float(np.mean(sample)) if len(sample) else None,
            peak_std=float(np.std(sample)) if len(sample) else None,
            moment_interval=moment_interval,
            interpretation='Local excess over equal-density adjacent sidebands; not a test of nonphysical origin'))
    return results


def conditional_response(truth_slope, truth_p, truth_x, truth_y, truth_entries,
                         reco_slope, reco_truth_index, reco_valid, reco_entries, split, edges, seed=50520,
                         bootstrap_repeats=100, windows=None, truth_selected=None,
                         phase_period=None):
    """Learn empirical residuals/multiplicity on discovery and predict held-out truth.

    Cells use predeclared p and transverse-position bins. The validation reco
    values enter only the final comparison, never the response or its weights.
    """
    truth_ok = np.isfinite(truth_slope) & np.isfinite(truth_p) & np.isfinite(truth_x) & np.isfinite(truth_y)
    if truth_selected is not None:
        truth_ok &= truth_selected
    momentum = np.searchsorted([50.,100.,250.,500.,1000.],truth_p)
    xbin = np.searchsorted([-30.,0.,30.],truth_x)
    ybin = np.searchsorted([-30.,0.,30.],truth_y)
    sbin = np.searchsorted([-.004,-.002,0.,.002,.004],truth_slope)
    cell = (momentum*16+xbin*4+ybin)*6+sbin
    if phase_period is not None:
        phase=np.floor(np.mod(np.where(np.isfinite(truth_slope),truth_slope,0.),phase_period)/phase_period*8).astype(int)
        cell=(((momentum*4+(truth_x>=0)*2+(truth_y>=0))*6+sbin)*8+phase)
    train_truth = truth_ok & (truth_entries<split)
    validation_truth = truth_ok & (truth_entries>=split)
    train_reco = reco_valid & np.isfinite(reco_slope) & (reco_entries<split)
    train_reco &= truth_ok[reco_truth_index]
    val_reco = reco_valid & np.isfinite(reco_slope) & (reco_entries>=split)
    val_reco &= truth_ok[reco_truth_index]
    rng = np.random.default_rng(seed)
    prediction = np.zeros(len(edges)-1)
    prediction2 = np.zeros(len(edges)-1)
    event_count = int(max(np.max(truth_entries,initial=-1),np.max(reco_entries,initial=-1))+1)
    # Poisson event bootstrap uses the same weight for all truths and tracks
    # from one event, including training efficiency and held-out observations.
    bootstrap_rng = np.random.default_rng(seed+1)
    event_weights = bootstrap_rng.poisson(1.,size=(bootstrap_repeats,event_count)).astype(np.uint8)
    predicted_bootstrap = np.zeros((bootstrap_repeats,len(edges)-1))
    observed_bootstrap = np.zeros_like(predicted_bootstrap)
    expected_total,underflow,overflow = 0.,0.,0.
    windows = windows or []
    predicted_moments = np.zeros((len(windows),3))
    predicted_bootstrap_moments = np.zeros((bootstrap_repeats,len(windows),3))
    observed_moments = np.zeros_like(predicted_moments)
    observed_bootstrap_moments = np.zeros_like(predicted_bootstrap_moments)
    def window_moments(values, weights):
        values = np.asarray(values).ravel()
        weights = np.broadcast_to(weights,values.shape)
        result = np.zeros((len(windows),3))
        for i,w in enumerate(windows):
            keep = np.abs(values-w['center'])<w['half_width']
            result[i] = [np.sum(weights[keep]),np.sum(weights[keep]*values[keep]),np.sum(weights[keep]*values[keep]**2)]
        return result
    zero_response_truth,zero_response_upper_sum = 0,0.
    cells, uncovered,covered_cells = [], 0,[]
    for c in np.unique(cell[validation_truth]):
        t = np.flatnonzero(train_truth & (cell==c))
        v = np.flatnonzero(validation_truth & (cell==c))
        r = np.flatnonzero(train_reco & (cell[reco_truth_index]==c))
        if len(t)>=30 and len(r)==0:
            covered_cells.append(c)
            zero_response_truth += len(v)
            zero_response_upper_sum += len(v)*(-np.log(.05))/len(t)
            cells.append(dict(cell=int(c),training_truth=len(t),training_reco=0,validation_truth=len(v),
                              covered=True,multiplicity=0.,zero_observed_training_response=True))
            continue
        if len(t)<30 or len(r)<20:
            uncovered += len(v)
            cells.append(dict(cell=int(c),training_truth=len(t),training_reco=len(r),validation_truth=len(v),covered=False))
            continue
        covered_cells.append(c)
        residuals = reco_slope[r]-truth_slope[reco_truth_index[r]]
        multiplicity = len(r)/len(t)
        # Eight independently drawn responses per validation truth; retain
        # under/overflow separately through the exported expected total.
        samples = truth_slope[v,None]+rng.choice(residuals,size=(len(v),8))
        weight = multiplicity/8
        prediction += np.histogram(samples,edges)[0]*weight
        prediction2 += np.histogram(samples,edges)[0]*weight*weight
        predicted_moments += window_moments(samples,weight)
        expected_total += len(v)*multiplicity
        underflow += np.count_nonzero(samples<edges[0])*weight
        overflow += np.count_nonzero(samples>=edges[-1])*weight
        for b,weights in enumerate(event_weights):
            truth_weight = weights[truth_entries[t]].sum()
            residual_weights = weights[reco_entries[r]].astype(float)
            if not truth_weight or not residual_weights.sum():
                continue
            multiplicity_b = residual_weights.sum()/truth_weight
            sampled = truth_slope[v,None]+bootstrap_rng.choice(residuals,size=(len(v),4),p=residual_weights/residual_weights.sum())
            sampled_weights = np.repeat(weights[truth_entries[v]]*multiplicity_b/4,4)
            predicted_bootstrap[b] += np.histogram(sampled.ravel(),edges,weights=sampled_weights)[0]
            predicted_bootstrap_moments[b] += window_moments(sampled,sampled_weights)
        cells.append(dict(cell=int(c),training_truth=len(t),training_reco=len(r),
                          validation_truth=len(v),covered=True,multiplicity=multiplicity,
                          residual_mean=float(np.mean(residuals)),residual_std=float(np.std(residuals))))
    observed = np.histogram(reco_slope[val_reco],edges)[0]
    uncovered_reco=val_reco & ~np.isin(cell[reco_truth_index],covered_cells)
    covered_reco=val_reco & ~uncovered_reco
    observed_covered=np.histogram(reco_slope[covered_reco],edges)[0]
    observed_covered_bootstrap=np.zeros_like(observed_bootstrap)
    covered_moments=window_moments(reco_slope[covered_reco],1.)
    covered_bootstrap_moments=np.zeros_like(observed_bootstrap_moments)
    observed_moments = window_moments(reco_slope[val_reco],1.)
    for b,weights in enumerate(event_weights):
        observed_bootstrap[b] = np.histogram(reco_slope[val_reco],edges,weights=weights[reco_entries[val_reco]].astype(float))[0]
        observed_bootstrap_moments[b] = window_moments(reco_slope[val_reco],weights[reco_entries[val_reco]].astype(float))
        observed_covered_bootstrap[b]=np.histogram(reco_slope[covered_reco],edges,weights=weights[reco_entries[covered_reco]].astype(float))[0]
        covered_bootstrap_moments[b]=window_moments(reco_slope[covered_reco],weights[reco_entries[covered_reco]].astype(float))
    a,b = prediction/max(1.,prediction.sum()),observed/max(1,observed.sum())
    midpoint = .5*(a+b)
    def divergence(p):
        good = p>0
        return np.sum(p[good]*np.log(p[good]/midpoint[good]))
    pa = predicted_bootstrap/np.maximum(1.,predicted_bootstrap.sum(axis=1))[:,None]
    pb = observed_bootstrap/np.maximum(1.,observed_bootstrap.sum(axis=1))[:,None]
    tv_bootstrap = .5*np.abs(pa-pb).sum(axis=1)
    centers = .5*(edges[1:]+edges[:-1])
    peak_comparison = []
    covered_peak_comparison=[]
    for i,window in enumerate(windows):
        predicted_n,observed_n = predicted_moments[i,0],observed_moments[i,0]
        def moments(values):
            n = values[...,0]
            mean = values[...,1]/np.maximum(1.,n)
            std = np.sqrt(np.maximum(0,values[...,2]/np.maximum(1.,n)-mean**2))
            return mean,std
        predicted_m,observed_m = moments(predicted_moments[i]),moments(observed_moments[i])
        boot_pred = moments(predicted_bootstrap_moments[:,i])
        boot_obs = moments(observed_bootstrap_moments[:,i])
        peak_comparison.append(dict(window=window,predicted_entries=float(predicted_n),observed_entries=int(observed_n),
            observed_minus_predicted_entries=float(observed_n-predicted_n),
            difference_entries_interval_95=np.quantile(observed_bootstrap_moments[:,i,0]-predicted_bootstrap_moments[:,i,0],[.025,.975]).tolist(),
            predicted_mean=float(predicted_m[0]),observed_mean=float(observed_m[0]),
            predicted_std=float(predicted_m[1]),observed_std=float(observed_m[1]),
            difference_mean_interval_95=np.quantile(boot_obs[0]-boot_pred[0],[.025,.975]).tolist(),
            difference_std_interval_95=np.quantile(boot_obs[1]-boot_pred[1],[.025,.975]).tolist(),
            moment_definition='Unbinned conditional moments inside frozen discovery windows'))
        covered_n=covered_moments[i,0]
        covered_m= moments(covered_moments[i])
        covered_boot=moments(covered_bootstrap_moments[:,i])
        covered_peak_comparison.append(dict(window=window,predicted_entries=float(predicted_n),
            observed_entries=int(covered_n),uncovered_observed_entries=int(observed_n-covered_n),
            observed_minus_predicted_entries=float(covered_n-predicted_n),
            difference_entries_interval_95=np.quantile(covered_bootstrap_moments[:,i,0]-predicted_bootstrap_moments[:,i,0],[.025,.975]).tolist(),
            predicted_mean=float(predicted_m[0]),observed_mean=float(covered_m[0]),
            predicted_std=float(predicted_m[1]),observed_std=float(covered_m[1]),
            difference_mean_interval_95=np.quantile(covered_boot[0]-boot_pred[0],[.025,.975]).tolist(),
            difference_std_interval_95=np.quantile(covered_boot[1]-boot_pred[1],[.025,.975]).tolist(),
            cohort='Only validation truths in covered training cells; same coverage for prediction and observation'))
    covered_a=prediction/max(1.,prediction.sum())
    covered_b=observed_covered/max(1,observed_covered.sum())
    covered_pa=predicted_bootstrap/np.maximum(1.,predicted_bootstrap.sum(axis=1))[:,None]
    covered_pb=observed_covered_bootstrap/np.maximum(1.,observed_covered_bootstrap.sum(axis=1))[:,None]
    return dict(predicted_counts=prediction.tolist(),mc_sumw2=prediction2.tolist(),
        predicted_counts_interval_95=np.quantile(predicted_bootstrap,[.025,.975],axis=0).tolist(),
        observed_counts_interval_95=np.quantile(observed_bootstrap,[.025,.975],axis=0).tolist(),
        observed_minus_predicted_interval_95=np.quantile(observed_bootstrap-predicted_bootstrap,[.025,.975],axis=0).tolist(),
        normalized_shape_total_variation=float(.5*np.abs(a-b).sum()),
        normalized_shape_total_variation_interval_95=np.quantile(tv_bootstrap,[.025,.975]).tolist(),
        prediction_expected_total=expected_total,prediction_underflow=underflow,prediction_overflow=overflow,
        peak_comparison=peak_comparison,
        covered_peak_comparison=covered_peak_comparison,
        covered_observed_counts=observed_covered.tolist(),
        covered_validation_reco=int(covered_reco.sum()),
        covered_observed_minus_predicted_interval_95=np.quantile(observed_covered_bootstrap-predicted_bootstrap,[.025,.975],axis=0).tolist(),
        covered_normalized_shape_total_variation=float(.5*np.abs(covered_a-covered_b).sum()),
        covered_normalized_shape_total_variation_interval_95=np.quantile(.5*np.abs(covered_pa-covered_pb).sum(axis=1),[.025,.975]).tolist(),
        bootstrap=dict(repeats=bootstrap_repeats,unit='source event',method='Poisson event bootstrap; four response draws per validation truth per replicate',
                       uncertainty='Training response/multiplicity, validation truth and observed tracks; includes finite response Monte Carlo variance'),
        observed_counts=observed.tolist(),edges=edges.tolist(),cells=cells,
        validation_truth=int(validation_truth.sum()),uncovered_validation_truth=uncovered,
        zero_training_response_validation_truth=zero_response_truth,
        zero_response_sum_per_cell_95_upper_expected_tracks=float(zero_response_upper_sum),
        uncovered_validation_reco=int(uncovered_reco.sum()),
        uncovered_validation_reco_counts=np.histogram(reco_slope[uncovered_reco],edges)[0].tolist(),
        discovery_reco=int(train_reco.sum()),validation_reco=int(val_reco.sum()),
        maximum_cdf_difference=float(np.max(np.abs(np.cumsum(a)-np.cumsum(b)),initial=0)),
        jensen_shannon_distance=float(np.sqrt(.5*(divergence(a)+divergence(b)))),
        model='Empirical residual and track multiplicity conditional on production p, same-plane truth slope, and true sensor x/y',
        residual_reference='Calypso-transported truth direction at the same plane; model reference',
        seed=seed,validation_reco_used_for_training=False,
        conditioning=dict(momentum_edges_GeV=[50,100,250,500,1000],
            coarse_slope_edges=[-.004,-.002,0.,.002,.004],
            transverse_position_edges_mm=[0.] if phase_period is not None else [-30.,0.,30.],
            phase_period=phase_period,phase_bins=8 if phase_period is not None else None))
