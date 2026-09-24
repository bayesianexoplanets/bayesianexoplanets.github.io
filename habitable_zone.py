"""Kopparapu et al. (2014) habitable-zone edges and the instellation of a planet, shared by the HZ figure
(make_habitable_zone.py) and the new-candidate table (build_new_candidates.py)."""
import numpy as np

# HZ-flux polynomial; i: 0 = Recent Venus, 1 = Runaway Greenhouse, 2 = Maximum Greenhouse, 3 = Early Mars.
SEFFSUN = [1.776, 1.107, 0.356, 0.320]
A = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5]
B = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9]
C = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12]
D = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16]


def hz_edge(teff, i):
    """Effective flux S_eff of HZ edge i at stellar temperature teff [K]."""
    Ts = teff - 5780.0
    return SEFFSUN[i] + A[i] * Ts + B[i] * Ts**2 + C[i] * Ts**3 + D[i] * Ts**4


def insolation(rstar, logg, teff, period):
    """Bolometric flux in Earth units (S/S_earth)."""
    mstar = 10.0**logg * rstar**2.0 / 10.0**4.437
    semia = mstar ** (1.0 / 3.0) * (period / 365.25) ** (2.0 / 3.0)
    lum = rstar**2.0 * (teff / 5778.0) ** 4.0
    return lum / semia**2.0


def zone(rstar, logg, teff, period):
    """'conservative' (runaway to maximum greenhouse), 'optimistic' (recent Venus to early Mars, outside the
    conservative zone), 'no' outside both, NaN without the stellar parameters."""
    s = insolation(rstar, logg, teff, period)
    if not np.isfinite(s) or not np.isfinite(teff):
        return np.nan
    if hz_edge(teff, 2) <= s <= hz_edge(teff, 1):
        return "conservative"
    if hz_edge(teff, 3) <= s <= hz_edge(teff, 0):
        return "optimistic"
    return "no"
