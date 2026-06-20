"""Singh-Maddala p-value helper.

The null/NST significance is the posterior-mean survival function of a per-star Singh-Maddala
(Burr XII) fit (informed log-normal priors; see overview.md + pvalue/fit_all.py). That posterior-mean
SF is a mixture, so it is stored per star as log10(SF) tabulated on a shared SNR grid (the
'sm_sf_grid' column). The p-value at any SNR is an interpolation of that grid:

    log10(p) = interp( log10(SNR);  log10(XGRID),  sf_grid ) .
"""
import numpy as np

# dense to 20 (nulls + steep SF transitions), coarse tail to 2500. MUST match pvalue/fit_all.py + index.html.
XGRID = np.unique(np.concatenate([np.geomspace(1.5, 20.0, 65), np.geomspace(20.0, 2500.0, 16)]))
LXG = np.log10(XGRID)


def log10p_from_grid(sf_grid_str, snr):
    """log10(p-value) at SNR from a star's pipe-separated log10-SF grid string. NaN if unusable."""
    try:
        g = np.array([float(v) for v in str(sf_grid_str).split("|")], dtype=float)
    except (ValueError, AttributeError):
        return np.nan
    if g.size != XGRID.size or not np.isfinite(snr) or snr <= 0:
        return np.nan
    return float(np.interp(np.log10(snr), LXG, g))   # grid is monotone; np.interp clamps at the ends
