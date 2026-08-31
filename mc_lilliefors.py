import numpy as np
from statsmodels.stats.diagnostic import lilliefors


def sample_truncated_exponential(beta, delta_m, size, rng=None):
    """
    Sample from the truncated exponential defined on the support [0, ∆m)
    """
    if rng is None:
        rng = np.random.default_rng()
    U = rng.random(size)
    return - (1/beta) * np.log(1 - U * (1 - np.exp(-beta * delta_m)))


def dither_magnitudes(mags, delta_m, b_value, rng=None):
    """
    Apply truncated-exponential dithering to binned magnitudes
    ----------
    mags : array of binned magnitudes
    delta_m : bin width 
    b_value : GR b-value
    -------
    Returns continuous magnitudes after dithering
    """
    MinMag = min(mags)
    beta = b_value * np.log(10)
    noise = sample_truncated_exponential(beta, delta_m, size=len(mags), rng=rng)
    return mags + noise - MinMag


def lilliefors_pval(mags):
    """
    Compute Lilliefors p-value for exponentiality
    """
    _, p_val = lilliefors(mags, dist='exp', pvalmethod='table')
    return p_val


MIN_EVENTS = 50
def estimate_mc_lilliefors(
    magnitudes,
    delta_m=0.1,
    b_value=1.0,
    alpha=0.1,
    n_dithers=50,
    min_events=MIN_EVENTS,
    rng=None,
):
    """
    Estimate magnitude of completeness Mc using the Lilliefors test with truncated-exponential dithering
    ----------
    magnitudes : array of binned magnitudes
    delta_m : bin width
    b_value : GR b-value
    alpha : significance level for the Lilliefors test
    n_dithers : number of dithering realizations (to calculate average p-value)
    min_events : minimum events required above Mc
    -------
    Returns:
    mc : estimated magnitude of completeness
    pvals : (dictionary) Mc candidate -> mean Lilliefors p-value
    """

    if not np.isfinite(delta_m) or delta_m <= 0:
        raise ValueError("delta_m must be positive and finite.")
    if not np.isfinite(b_value) or b_value <= 0:
        raise ValueError("b_value must be positive and finite.")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1.")
    if n_dithers < 1:
        raise ValueError("n_dithers must be at least 1.")
    if min_events < 1:
        raise ValueError("min_events must be at least 1.")

    mags = np.asarray(magnitudes, dtype=float)
    mags = mags[np.isfinite(mags)]
    if mags.size == 0:
        raise ValueError("No finite magnitudes were provided.")
    unique_bins = np.sort(np.unique(mags))
    if rng is None:
        rng = np.random.default_rng()

    pval_by_mc = {}

    for mc in unique_bins:
        subset = mags[mags >= mc]
        if len(subset) < min_events:
            continue

        pvals = []
        for _ in range(n_dithers):
            dith_mags = dither_magnitudes(subset, delta_m, b_value, rng=rng)
            pvals.append(lilliefors_pval(dith_mags))

        pval_by_mc[mc] = np.mean(pvals)

    mc_candidates = [m for m, p in pval_by_mc.items() if p > alpha]

    if len(mc_candidates) == 0:
        return None, pval_by_mc 
    
    else:
        # Mc = smallest magnitude bin where the null is not rejected
        mc_est = mc_candidates[0]
        return mc_est, pval_by_mc
