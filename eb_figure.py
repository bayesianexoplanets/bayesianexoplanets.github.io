"""eb_figure.py -- folded figures of the website's eclipsing binaries with a model curve, on the all-sky archive.

Every row of `tess_ebs.csv` on a star of the archive gets one figure: the primary and the secondary eclipse on top, the
whole orbit below (a contact binary: the whole orbit only). The dashed curve is, by row:

- posterior : rows moved here by the false-positive analysis. The pointwise median of the model curves of the samples
              of the leading binary scenario of the Pentaceratops posterior, with the 16-84 % band; only a constant
              baseline and a time shift are fitted (and one dilution factor when the binary is a resolved nearby
              star). The fold uses the SPOC two-minute sectors only, the data that posterior was fitted to: QLP
              sectors change the depth of crowded stars by up to a quarter.
- fit       : every other row. The same two-sphere eclipse model fitted to the fold, started from a two-Gaussian fit;
              the two-Gaussian curve itself for contact binaries and where it fits better. Fitted parameters are
              degenerate and are not published.
- none      : the data only: rows reviewed as "no signal" whose fold shows no eclipse either, and rows whose fitted
              curve does not describe the fold (`eb_routes.tsv`, written from the records of a first pass).

The model is TRICERATOPS's eclipsing-binary light curve (a limb-darkened eclipse of the host by a sphere, diluted by the
other light) with the secondary eclipse added; design, checks and sources: logs/DESIGN.md 2026-11-19 and
`allsky_search/eb_theory_curve/PROPOSAL.md`. The light curve is detrended in two passes so that the GP never bends the
eclipses: a first GP on the flux minus the whole-orbit template, then a GP on the flux with the measured eclipse windows
masked. Units are the pipeline's (the light curve over its clipped standard deviation); the error bars carry the
measured scatter around the curve, because for a binary that standard deviation is set by the eclipses, not the noise.

    python eb_figure.py list                           the rows, by route
    python eb_figure.py draw [workers] [i/n] [TIC ...]  figures to scratch, resumable; i/n takes every n-th row from i
    python eb_figure.py summary                        what was drawn, how well the curves fit
    python eb_figure.py install [--apply]              WebP files into plots_eb/, results into eb_figures.tsv (dry run
                                                       by default); then run build_verified_catalogs.py

Runs in the pipeline's environment (`cosmology`, Python 3.8).
"""
import os
import sys
import time
import warnings

sys.dont_write_bytecode = True
for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from constants import scratch, dt
import site_figures

ARCHIVE = scratch + 'allsky/archive'
OUT = scratch + 'EB_figures_20261119/'
EB_TABLE = os.path.join(HERE, 'tess_ebs.csv')
RESULTS = os.path.join(HERE, 'eb_figures.tsv')         # read by build_verified_catalogs.py
ROUTES = os.path.join(HERE, 'eb_routes.tsv')           # Name, route: rows whose route is set by hand or by a first pass
RAW_CATALOG = scratch + 'allsky_search/known_objects/raw_20261114/tess_ebs.csv'
POSTERIOR = scratch + 'TRICERATOPS/pentaceratops_fpp_20261117/posterior_false_positive.npz'

EXPOSURE_LONG = 0.5 / 24.        # 30-min FFI integration [d]
SUPERSAMPLE = 15                 # batman sub-samples per exposure
ELSEWHERE = ('STP', 'BTP', 'NTP')                      # planet-on-another-star scenarios: not binaries
COMPANION_HOST = ('SEB', 'SEBx2P', 'BEB', 'BEBx2P')    # the eclipsing pair orbits the unresolved companion
NEARBY_HOST = ('NEB', 'NEBx2P')                        # the pair is a resolved nearby star: aperture dilution unknown
MAX_MASKED = 0.45                # largest share of valid samples a GP mask may remove
CONTACT_WIDTH = 0.2              # eclipse width over period above which only the whole orbit is drawn
CONTACT_MORPH = 0.7              # TESS-EB morphology parameter above which a binary counts as contact
N_POSTERIOR_CURVES = 300         # posterior samples drawn per figure
PHYSICAL_MARGIN = 1.2            # the two-sphere fit is kept while its chi^2 is below this times the two-Gaussian's
NO_SIGNAL = 'no signal'          # ending of the `Checked` label of rows where a review found nothing
SPOC_ORIGINS = (0,)              # TESS/shard_archive.ORIGIN_CODES: the SPOC two-minute light curves
MAX_SHIFT = 0.25                 # largest fitted time shift of a posterior curve, in listed eclipse widths


# ----------------------------------------------------------------------------------------------------------------------
# Inputs

class EBRows:
    """The website's EB table with what each figure needs: its file, its route, catalog shape columns, posterior samples."""

    def __init__(self):
        table = pd.read_csv(EB_TABLE, low_memory=False)
        table['index'] = table.groupby('TIC').cumcount()
        table['figure'] = 'plots_eb/' + table['TIC'].astype(str) + '/' + table['index'].astype(str) + '.webp'
        table['route'] = np.where(table['Checked'].str.endswith(NO_SIGNAL), 'none',
                                  np.where(table['Origin'] == 'false-positive analysis', 'posterior', 'fit'))
        self.table = table
        self.catalog = self.catalog_shapes()
        self.posterior = PosteriorFile()
        missing = (table['route'] == 'posterior') & ~table['Name'].isin(self.posterior.index)
        table.loc[missing, 'route'] = 'fit'               # no samples for this row: fitted like a catalog binary
        if os.path.exists(ROUTES):
            forced = pd.read_csv(ROUTES, sep='\t').set_index('Name')['route']
            table['route'] = table['Name'].map(forced).fillna(table['route'])


    @staticmethod
    def catalog_shapes():
        """Morphology and two-Gaussian eclipse widths of the TESS-EB catalog (Prsa et al. 2022), by catalog entry."""
        raw = pd.read_csv(RAW_CATALOG, comment='#', sep='\t', skiprows=[1, 2], skipinitialspace=True)
        for column in raw.columns:
            raw[column] = pd.to_numeric(raw[column], errors='coerce')
        raw = raw.dropna(subset=['TIC'])
        raw['entry'] = ['%010d.%d' % (tic, k + 1) for tic, k in zip(raw['TIC'].astype(int), raw.groupby('TIC').cumcount())]
        return raw.set_index('entry')[['Morph', 'Wp-2g', 'Ws-2g']]


    def tasks(self, only=None, part=(0, 1)):
        """One picklable task per row on an archive star."""
        rows = self.table[self.table['in_archive'] == True]
        rows = rows[rows['TIC'].isin(only)] if only else rows.iloc[part[0]::part[1]]
        out = []
        for r in rows.itertuples():
            entry = str(r.Catalog) if isinstance(r.Catalog, str) else ''
            shape = self.catalog.loc[entry].to_dict() if entry in self.catalog.index else {}
            out.append(dict(name=r.Name, tic=int(r.TIC), figure=r.figure, route=r.route, period=float(r.Period), epoch=float(r.Epoch),
                            duration=float(r.Duration), secondary_phase=float(r.Secondary_phase),
                            secondary_duration=float(r.Secondary_duration), morph=shape.get('Morph', np.nan),
                            width_1=shape.get('Wp-2g', np.nan), width_2=shape.get('Ws-2g', np.nan),
                            posterior=self.posterior.leading(r.Name) if r.route == 'posterior' else None))
        return out


class PosteriorFile:
    """The Pentaceratops false-positive posterior (its README): per candidate, samples of every scenario."""

    def __init__(self, path=POSTERIOR):
        with np.load(path, allow_pickle=False) as archive:
            self.names = archive['candidate'].astype(str)
            self.start, self.stop = archive['start'], archive['stop']
            self.columns = [str(c) for c in archive['columns']]
            self.scenario_names = np.array([str(s) for s in archive['scenario_names']])
            self.samples, self.scenario = archive['samples'], archive['scenario']
        self.index = {name: i for i, name in enumerate(self.names)}


    def leading(self, name, seed=1):
        """The most populated binary scenario of a candidate: its name, its share of all samples, N_POSTERIOR_CURVES of
        its samples and the mean epoch offset of the binary samples (the table epoch is the catalogue epoch plus it)."""
        i = self.index[name]
        frame = pd.DataFrame(self.samples[self.start[i]:self.stop[i]].astype(float), columns=self.columns)
        frame['scenario'] = self.scenario_names[self.scenario[self.start[i]:self.stop[i]]]
        binary = frame[~frame['scenario'].isin(ELSEWHERE)]
        counts = binary['scenario'].value_counts()
        samples = binary[binary['scenario'] == counts.index[0]]
        take = np.random.default_rng(seed).choice(len(samples), size=min(N_POSTERIOR_CURVES, len(samples)), replace=False)
        return dict(scenario=counts.index[0], share=counts.iloc[0] / len(frame), n_samples=len(samples),
                    mean_offset=float(binary['epoch_offset_min'].mean()),
                    samples=samples.iloc[np.sort(take)].drop(columns='scenario').reset_index(drop=True))


# ----------------------------------------------------------------------------------------------------------------------
# Light curve

class DetrendedLightCurve:
    """A star's archive light curve in pipeline units with the pipeline GP removed, the eclipses left intact.

    Pass 1 fits the GP to the flux minus the whole-orbit template (outlier filter and a first detrending that never
    sees the eclipses); the two-Gaussian model of that fold measures the eclipse centres and widths w; pass 2 refits
    the GP to the flux with windows of +-0.75 w around both eclipses masked, so it follows the out-of-eclipse
    variability and interpolates only across the measured eclipses. Windows over MAX_MASKED of the data (contact
    binaries, wide shallow eclipses): pass 1 is kept.

    Parameters
    ----------
    tic : int
    period_orb, epoch : float
        Orbital period [d] and time of the primary eclipse [BTJD].
    hints : list of (float, float)
        Starting (offset from the epoch [d], half-width [d]) of the primary and the secondary eclipse.
    origins : tuple of int, optional
        Block origins the fold keeps (every block is detrended); all blocks when None or when the star has none of them.
    """

    def __init__(self, tic, period_orb, epoch, hints, origins=None):
        from pipeline import fold_plot
        from TESS.load_tess import load_stellar_data, load_origin, StarInfo_tess
        self.tic, self.period_orb, self.epoch = tic, period_orb, epoch
        stellar = load_stellar_data(tic)
        self.star = StarInfo_tess(*stellar[1:])
        self.sigma = float(self.star.sigma_star_flux_units)
        Flux, invVar, Indices, upsampled = fold_plot.get_data(tic)
        origin = load_origin(tic).origin
        btjd = float(stellar[0]) + Indices * dt
        valid = invVar > 0.5
        folded = valid & np.isin(np.asarray(origin), origins)[:, None] if origins is not None else valid
        folded = folded if folded.any() else valid

        template = np.where(valid, self.template(btjd, Flux, invVar, folded), 0.)
        gflux, fgp, _ = fold_plot.iteration(upsampled, Flux - template, invVar, Indices, origin)
        filtered = gflux + template                          # outlier-filtered flux, eclipses intact
        data = filtered - fgp
        self.windows = self.measure_windows(btjd, data, invVar, folded, hints)
        masked = self.window_mask(btjd, self.windows, valid, min_scale=2. / 3.)   # down to +-0.5 w
        self.masked_share = float(masked[valid].mean())
        if self.masked_share <= MAX_MASKED:                  # detached: GP refitted with the measured windows masked
            _, fgp, _ = fold_plot.iteration(upsampled, Flux, np.where(masked, 0., invVar), Indices, origin)
            data = filtered - fgp
        self.time, self.flux, self.weight = btjd[folded], data[folded], invVar[folded]
        self.exposure = EXPOSURE_LONG if np.any(upsampled) else dt
        self.n_blocks, self.n_folded = int(Flux.shape[0]), int(folded.any(axis=1).sum())


    def template(self, btjd, flux, invVar, valid):
        """Weighted phase-binned mean of the flux at the orbital period, one bin per pipeline sample, empty bins
        interpolated periodically; returned at the samples."""
        n_bins = max(int(round(self.period_orb / dt)), 8)
        phase = np.mod(btjd - self.epoch, self.period_orb) / self.period_orb
        bins = np.minimum((phase * n_bins).astype(int), n_bins - 1)
        weight = np.where(valid, invVar, 0.)
        total = np.bincount(bins[valid], weights=weight[valid], minlength=n_bins)
        summed = np.bincount(bins[valid], weights=(flux * weight)[valid], minlength=n_bins)
        filled = total > 0
        centres = (np.arange(n_bins) + 0.5) / n_bins
        mean = np.interp(centres, centres[filled], summed[filled] / total[filled], period=1.)
        return mean[bins]


    def measure_windows(self, btjd, flux, invVar, valid, hints):
        """Eclipse windows on the fold of the first pass: the two-Gaussian model's centres and widths w = 5.6 sigma
        (Mowlavi et al. 2017, Eqs. 2 and 4), half-window 0.75 w."""
        P = self.period_orb
        weight = np.where(valid, invVar, 0.)
        x = np.mod(btjd[valid] - self.epoch + 0.25 * P, P) - 0.25 * P
        n_bins = max(int(round(P / dt)), 8)
        bins = np.minimum(((x + 0.25 * P) / P * n_bins).astype(int), n_bins - 1)
        total = np.bincount(bins, weights=weight[valid], minlength=n_bins)
        summed = np.bincount(bins, weights=(flux * weight)[valid], minlength=n_bins)
        filled = total > 0.3
        centres = -0.25 * P + (np.arange(n_bins) + 0.5) * P / n_bins
        offset_2 = np.mod(hints[1][0] + 0.25 * P, P) - 0.25 * P
        width = 2. * hints[0][1] / 1.5 / P
        fit = TwoGaussianFit(centres[filled], summed[filled] / total[filled], 1. / np.sqrt(total[filled]), P,
                             mu_2=offset_2 / P, width_1=width, width_2=width)
        fit.run()
        v = fit.minuit.values
        widths = [min(5.6 * v['s_1'], 0.4) * P, min(5.6 * v['s_2'], 0.4) * P]
        if fit.coef[2] < 3. * fit.coef_err[2]:               # no significant secondary: the primary's width at the hint
            widths[1], mu_2 = widths[0], offset_2 / P
        else:                                                # a shallow secondary must not take the orbit
            widths[1], mu_2 = min(widths[1], 2. * widths[0]), v['mu_2']
        return [(v['mu_1'] * P, 0.75 * widths[0]), (mu_2 * P, 0.75 * widths[1])]


    def window_mask(self, btjd, windows, valid, min_scale=0.05):
        """Samples inside any (centre offset from the epoch [d], half-width [d]) window; all half-widths shrink
        together (factor 0.8 per step, not below min_scale) until at most MAX_MASKED of the valid samples is masked."""
        scale = 1.
        while True:
            masked = np.zeros(btjd.shape, dtype=bool)
            for offset, half in windows:
                x = np.mod(btjd - self.epoch - offset + 0.5 * self.period_orb, self.period_orb) - 0.5 * self.period_orb
                masked |= np.abs(x) < scale * half
            if masked[valid].mean() <= MAX_MASKED or scale * 0.8 < min_scale:
                return masked
            scale *= 0.8


    def phase_time(self):
        """Every sample's time from the primary eclipse, in [-P/4, 3P/4) [d]."""
        return np.mod(self.time - self.epoch + 0.25 * self.period_orb, self.period_orb) - 0.25 * self.period_orb


    def fold(self, bin_width):
        """Weighted fold at the orbital period from phase -0.25 over one period.

        Returns
        -------
        x : ndarray
            Bin centres, days from the primary eclipse.
        y, err : ndarray
            Weighted mean flux and 1 / sqrt(sum of weights), pipeline units.
        """
        P = self.period_orb
        x = self.phase_time()
        n_bins = max(int(np.ceil(P / bin_width)), 1)
        bins = np.floor((x + 0.25 * P) / bin_width).astype(int)
        keep = (bins >= 0) & (bins < n_bins)
        total = np.bincount(bins[keep], weights=self.weight[keep], minlength=n_bins)
        summed = np.bincount(bins[keep], weights=(self.flux * self.weight)[keep], minlength=n_bins)
        centres = -0.25 * P + (np.arange(n_bins) + 0.5) * bin_width
        filled = total > 0.3
        return centres[filled], summed[filled] / total[filled], 1. / np.sqrt(total[filled])


# ----------------------------------------------------------------------------------------------------------------------
# Models

class BinaryEclipseModel:
    """Two limb-darkened spheres on a Keplerian orbit, each eclipse a batman transit, light-weighted with any third
    light (the dilution of TRICERATOPS's simulate_EB_transit, with the secondary eclipse added).

    Parameters
    ----------
    exposure : float
        Integration time [d] over which batman averages the light curve.
    """

    def __init__(self, exposure):
        self.exposure = exposure


    def transit(self, t, t0, per, rp, a, inc, ecc, w, u):
        import batman
        params = batman.TransitParams()
        params.t0, params.per, params.rp, params.a, params.inc = t0, per, rp, a, inc
        params.ecc, params.w, params.u, params.limb_dark = ecc, w, [float(u[0]), float(u[1])], 'quadratic'
        if self.exposure > 0.:
            return batman.TransitModel(params, t, supersample_factor=SUPERSAMPLE, exp_time=self.exposure).light_curve(params)
        return batman.TransitModel(params, t).light_curve(params)


    @staticmethod
    def secondary_time(t0, per, rp, a, inc, ecc, w):
        """Time of the secondary eclipse (batman get_t_secondary), within one period after t0."""
        import batman
        params = batman.TransitParams()
        params.t0, params.per, params.rp, params.a, params.inc = t0, per, rp, a, inc
        params.ecc, params.w, params.u, params.limb_dark = ecc, w, [0.3, 0.2], 'quadratic'
        params.fp, params.t_secondary = 0.1, t0 + 0.5 * per
        t_sec = batman.TransitModel(params, np.array([t0])).get_t_secondary(params)
        return t0 + np.mod(t_sec - t0, per)


    def relative_flux(self, t, per, t0, k, a_r1, inc, ecc, w, light_2, u_1, u_2, light_other=0.):
        """Total flux over its uneclipsed value, (lambda_1 + L2 lambda_2 + L3) / (1 + L2 + L3): star 1 is eclipsed at
        the primary by star 2, star 2 at the secondary by star 1 (roles swapped: k -> 1/k, a/R1 -> a/R2,
        w -> w + 180 deg); luminosities in units of star 1's."""
        lam_1 = self.transit(t, t0, per, k, a_r1, inc, ecc, w, u_1)
        t_sec = self.secondary_time(t0, per, k, a_r1, inc, ecc, w)
        lam_2 = self.transit(t, t_sec, per, 1. / k, a_r1 / k, inc, ecc, w + 180., u_2)
        return (lam_1 + light_2 * lam_2 + light_other) / (1. + light_2 + light_other)


    def mean_flux(self, per, t0, n=3000, **args):
        """Orbit average of relative_flux on n uniform times, without exposure smearing (it leaves the mean unchanged)."""
        exposure, self.exposure = self.exposure, 0.
        mean = float(np.mean(self.relative_flux(t0 + np.arange(n) * per / n, per, t0, **args)))
        self.exposure = exposure
        return mean


def pipeline_units(relative, sigma, mean_relative):
    """Relative flux (1 out of eclipse) in the figure's units: each sector was divided by its mean flux, which the
    eclipses lower to mean_relative of the out-of-eclipse flux, then by sigma; a constant baseline is fitted apart."""
    return (relative / mean_relative - 1.) / sigma


def two_gaussian_basis(phase, mu_1, mu_2, s_1, s_2, phi_ell):
    """Columns of Mowlavi et al. (2017) Eq. 2 in flux: constant, the two Gaussian eclipses with their mirrors at
    m = -2..2 (minus sign: flux, not magnitude), and a cosine at half the period."""
    def gauss(mu, s):
        return sum(np.exp(-(phase - mu - m) ** 2 / (2. * s ** 2)) for m in range(-2, 3))
    return np.column_stack([np.ones_like(phase), -gauss(mu_1, s_1), -gauss(mu_2, s_2), 0.5 * np.cos(4. * np.pi * (phase - phi_ell))])


def linear_fit(y, err, basis):
    """Weighted linear least squares: coefficients, chi^2 and the coefficients' standard errors."""
    design = basis / err[:, None]
    coef = np.linalg.lstsq(design, y / err, rcond=None)[0]
    chi2 = float(np.sum((y / err - design @ coef) ** 2))
    return coef, chi2, np.sqrt(np.diag(np.linalg.pinv(design.T @ design)))


class TwoGaussianFit:
    """Mowlavi et al. (2017) Eq. 2 fitted to a fold: phases and widths of the eclipses by iminuit, the constant, both
    depths and the ellipsoidal amplitude by linear least squares at each step.

    Parameters
    ----------
    x, y, err : ndarray
        Fold: days from the primary, flux and error.
    period : float
    mu_2, width_1, width_2 : float
        Starting secondary phase and eclipse widths (phase units, w = 5.6 sigma).
    """

    def __init__(self, x, y, err, period, mu_2=0.5, width_1=0.05, width_2=0.05):
        self.period = period
        self.phase, self.y, self.err = x / period, y, err
        self.start = dict(mu_2=float(np.clip(mu_2, 0.2, 0.8)), s_1=max(width_1 / 5.6, 0.003), s_2=max(width_2 / 5.6, 0.003))


    def chi2(self, mu_1, mu_2, s_1, s_2, centre):
        phi_ell = mu_1 if centre == 1 else mu_2
        return linear_fit(self.y, self.err, two_gaussian_basis(self.phase, mu_1, mu_2, s_1, s_2, phi_ell))[1]


    def run(self):
        from iminuit import Minuit
        best = None
        for centre in (1, 2):
            m = Minuit(lambda mu_1, mu_2, s_1, s_2: self.chi2(mu_1, mu_2, s_1, s_2, centre), mu_1=0., **self.start)
            m.errordef = Minuit.LEAST_SQUARES
            m.limits['mu_1'], m.limits['mu_2'] = (-0.05, 0.05), (0.2, 0.8)
            m.limits['s_1'], m.limits['s_2'] = (0.001, 0.25), (0.001, 0.25)
            m.migrad()
            if best is None or m.fval < best[0].fval:
                best = (m, centre)
        self.minuit, self.centre = best
        self.coef, _, self.coef_err = linear_fit(self.y, self.err, self.basis(self.phase))
        return self.minuit


    def basis(self, phase):
        v = self.minuit.values
        return two_gaussian_basis(phase, v['mu_1'], v['mu_2'], v['s_1'], v['s_2'], v['mu_1'] if self.centre == 1 else v['mu_2'])


    def curve(self, x):
        return self.basis(x / self.period) @ self.coef


    def windows(self):
        """(centre [d], half-window 0.75 w [d]) of both eclipses, w = min(5.6 sigma, 0.4) (Mowlavi et al. 2017 Eq. 4)."""
        v = self.minuit.values
        return [(v['mu_1'] * self.period, 0.75 * min(5.6 * v['s_1'], 0.4) * self.period),
                (v['mu_2'] * self.period, 0.75 * min(5.6 * v['s_2'], 0.4) * self.period)]


class PhysicalFit:
    """BinaryEclipseModel fitted to the fold (iminuit), started from the two-Gaussian fit of the same fold. Free:
    shift, k = R2/R1, r_sum = (R1 + R2)/a, primary impact parameter b, light ratio L2/L1, e cos w, e sin w, third-light
    fraction, baseline. Fixed: period, limb darkening (the star's u1, u2 for both stars). The third light only lets
    shallow eclipses be reached; it is degenerate with the rest, not a measurement."""

    def __init__(self, star, exposure, x, y, err, period, windows):
        self.period = period
        self.u = (star.u1, star.u2)
        self.sigma = float(star.sigma_star_flux_units)
        self.model = BinaryEclipseModel(exposure)
        near = np.zeros(len(x), dtype=bool)                  # bins near either eclipse, plus every 8th bin elsewhere
        for centre, half in windows:
            near |= np.abs(np.mod(x - centre + 0.5 * period, period) - 0.5 * period) < 2. * half
        keep = near | (np.arange(len(x)) % 8 == 0)
        self.x, self.y, self.err = x[keep], y[keep], err[keep]


    @staticmethod
    def geometry(k, r_sum, b, ecosw, esinw):
        """(a/R1, inclination [deg], e, w [deg]) from the fitted shape parameters; b by Winn (2010) Eq. 7."""
        ecc = min(np.hypot(ecosw, esinw), 0.9)
        w = np.degrees(np.arctan2(esinw, ecosw))
        a_r1 = (1. + k) / r_sum
        cos_i = b * (1. + esinw) / (a_r1 * (1. - ecc ** 2))
        return a_r1, np.degrees(np.arccos(np.clip(cos_i, 0., 1.))), ecc, w


    def curve(self, x, shift, k, r_sum, b, light_2, ecosw, esinw, third, baseline):
        a_r1, inc, ecc, w = self.geometry(k, r_sum, b, ecosw, esinw)
        args = dict(per=self.period, t0=shift, k=k, a_r1=a_r1, inc=inc, ecc=ecc, w=w, light_2=light_2, u_1=self.u, u_2=self.u)
        relative = (1. - third) * self.model.relative_flux(x, **args) + third
        mean = (1. - third) * self.model.mean_flux(**args) + third
        return pipeline_units(relative, self.sigma, mean) + baseline


    def starts(self, gauss):
        """Starting points from the two-Gaussian fit: r_sum from the primary width, e cos w from the secondary phase
        (Winn 2010 Eq. 33), e sin w from the width ratio (Eq. 34), L2/L1 = (d2/d1) k^2, third light matching the depth."""
        v = gauss.minuit.values
        w_1, w_2 = min(5.6 * v['s_1'], 0.4), min(5.6 * v['s_2'], 0.4)
        d_1, d_2 = max(gauss.coef[1], 1e-6), max(gauss.coef[2], 0.)
        r_sum = float(np.clip(np.sin(np.pi * w_1), 0.01, 0.8))
        ecosw = float(np.clip(np.pi / 4. * (2. * (v['mu_2'] - v['mu_1']) - 1.), -0.5, 0.5))
        ratio = w_2 / w_1 if d_2 > 0 else 1.
        esinw = float(np.clip((ratio - 1.) / (ratio + 1.), -0.3, 0.3))
        depth = d_1 * self.sigma                             # relative depth of the primary
        starts = []
        for k in (0.3, 0.8):
            for b in (0.2, 0.6):
                light_2 = float(np.clip(d_2 / d_1 * k ** 2, 0.001, 5.))
                third = float(np.clip(1. - depth / (min(k ** 2, 1.) / (1. + light_2)), 0., 0.95))
                starts.append(dict(shift=v['mu_1'] * self.period, k=k, r_sum=r_sum, b=b, light_2=light_2, ecosw=ecosw,
                                   esinw=esinw, third=third, baseline=gauss.coef[0]))
        return starts


    def run(self, gauss, duration):
        from iminuit import Minuit
        from iminuit.cost import LeastSquares
        best = None
        for start in self.starts(gauss):
            m = Minuit(LeastSquares(self.x, self.y, self.err, self.curve), **start)
            m.limits['shift'] = (start['shift'] - 0.5 * duration, start['shift'] + 0.5 * duration)
            m.limits['k'], m.limits['r_sum'], m.limits['b'] = (0.05, 3.), (0.002, 0.9), (0., 2.)
            m.limits['light_2'], m.limits['ecosw'], m.limits['esinw'] = (0., 20.), (-0.6, 0.6), (-0.6, 0.6)
            m.limits['third'] = (0., 0.99)
            m.migrad()
            if best is None or m.fval < best.fval:
                best = m
        self.minuit = best
        return best


    def secondary_time(self):
        v = self.minuit.values
        a_r1, inc, ecc, w = self.geometry(v['k'], v['r_sum'], v['b'], v['ecosw'], v['esinw'])
        return self.model.secondary_time(v['shift'], self.period, v['k'], a_r1, inc, ecc, w)


class PosteriorCurve:
    """Model curves of the drawn samples of the leading binary scenario; every parameter from the posterior.

    Parameters
    ----------
    posterior : dict
        `PosteriorFile.leading` of the row.
    period, epoch : float
        The row's listed period and epoch.
    star : StarInfo_tess
        Supplies u1, u2 (used for both stars, as TRICERATOPS does) and sigma.
    exposure : float
    """
    GM_SUN = 2942.2062175044193          # G Msun / Rsun^3 [1 / day^2] (astropy constants GM_sun, R_sun)

    def __init__(self, posterior, period, epoch, star, exposure):
        self.scenario, self.share, self.samples = posterior['scenario'], posterior['share'], posterior['samples']
        self.t0 = epoch + (self.samples['epoch_offset_min'].values - posterior['mean_offset']) / 1440.
        self.period_orb = self.multiple(posterior, period) * period    # the table period; the file's float32 P_orb is coarser
        self.model = BinaryEclipseModel(exposure)
        self.u = (star.u1, star.u2)
        self.sigma = float(star.sigma_star_flux_units)
        self.nearby = self.scenario in NEARBY_HOST


    @staticmethod
    def multiple(posterior, period):
        """Orbital period over listed period: 2 in the x2P scenarios."""
        return int(round(np.median(posterior['samples']['P_orb']) / period))


    def parameters(self, i):
        s = self.samples.iloc[i]
        a_r1 = np.cbrt(self.GM_SUN * (s['M_s'] + s['M_EB']) * s['P_orb'] ** 2 / (4. * np.pi ** 2)) / s['R_s']
        f_eb = s['fluxratio_EB'] / (1. - s['fluxratio_EB'])
        f_comp = s['fluxratio_comp'] / (1. - s['fluxratio_comp'])
        if self.scenario in COMPANION_HOST:                  # host = the companion: weights (F_comp, F_EB, F_target = 1)
            light_2, light_other = f_eb / f_comp, 1. / f_comp
        else:                                                # host = target (or the nearby star, in its own units)
            light_2, light_other = f_eb, f_comp
        return dict(per=s['P_orb'], t0=self.t0[i], k=s['R_EB'] / s['R_s'], a_r1=a_r1, inc=s['inc'], ecc=s['ecc'], w=s['omega'],
                    light_2=light_2, light_other=light_other)


    def curves(self, x_grid, epoch):
        """Model curves of every drawn sample on x_grid (days from the epoch), pipeline units, shape (n, len)."""
        out = np.empty((len(self.samples), len(x_grid)))
        for i in range(len(self.samples)):
            p = self.parameters(i)
            relative = self.model.relative_flux(epoch + x_grid, u_1=self.u, u_2=self.u, **p)
            out[i] = pipeline_units(relative, self.sigma, self.model.mean_flux(u_1=self.u, u_2=self.u, **p))
        return out


    def secondary_offset(self):
        """Median secondary-eclipse time over the drawn samples, days after the epoch."""
        times = []
        for i in range(len(self.samples)):
            p = self.parameters(i)
            times.append(self.model.secondary_time(p['t0'], p['per'], p['k'], p['a_r1'], p['inc'], p['ecc'], p['w']) - self.t0.mean())
        return float(np.median(times))


# ----------------------------------------------------------------------------------------------------------------------
# Figure

class EBFigure:
    """One row's light curve, curve and figure.

    Parameters
    ----------
    task : dict
        A row as `EBRows.tasks` gives it.
    """

    def __init__(self, task):
        self.task = task
        self.record = dict(name=task['name'], TIC=task['tic'], figure=task['figure'], route=task['route'])
        self.curve, self.band, self.label, self.title = None, None, '', ''
        self.physical = True                                 # False: the two-Gaussian curve only (quick looks)


    def run(self, png):
        self.detrend()
        return self.finish(png)


    def detrend(self):
        """The row's light curve at its orbital period (twice the listed one for an x2P posterior)."""
        task = self.task
        os.environ['TESS_SHARD_ARCHIVE'] = ARCHIVE
        from TESS.load_tess import load_stellar_data, StarInfo_tess
        posterior = None
        self.period_orb = task['period']
        if task['route'] == 'posterior':
            posterior = PosteriorCurve(task['posterior'], task['period'], task['epoch'], StarInfo_tess(*load_stellar_data(task['tic'])[1:]), dt)
            self.period_orb = posterior.period_orb
        P, duration = self.period_orb, task['duration']

        hints = [(0., 1.5 * duration)]
        if np.isfinite(task['secondary_phase']):
            hints.append((task['secondary_phase'] * task['period'], 1.5 * np.nanmax([task['secondary_duration'], duration])))
        else:
            hints.append((posterior.secondary_offset() if posterior is not None else 0.5 * P, 1.5 * duration))
        start = time.time()
        self.lc = DetrendedLightCurve(task['tic'], P, task['epoch'], hints, SPOC_ORIGINS if posterior is not None else None)
        self.record.update(period_orb=P, n_blocks=self.lc.n_blocks, n_folded=self.lc.n_folded, masked_share=self.lc.masked_share,
                           detrend_s=time.time() - start)
        self.wide = duration / P > CONTACT_WIDTH or (np.isfinite(task['morph']) and task['morph'] >= CONTACT_MORPH)
        self.posterior = posterior


    def finish(self, png):
        """Curve, fit statistics and figure on the light curve `self.lc`."""
        task = self.task
        start = time.time()
        self.t_sec = (task['secondary_phase'] if np.isfinite(task['secondary_phase']) else 0.5) * task['period']
        if self.posterior is not None:
            self.from_posterior(self.posterior)
        elif task['route'] == 'fit':
            self.from_fit()
        self.place_secondary()
        self.record.update(model_s=time.time() - start, **self.goodness())
        self.draw(png)
        return self.record


    def from_posterior(self, posterior):
        """Median and 16-84 % band of the sample curves; a constant baseline and a time shift of at most MAX_SHIFT listed
        widths are fitted (nearby host: also one dilution factor)."""
        from scipy.optimize import minimize_scalar
        P, lc = self.period_orb, self.lc
        posterior.model.exposure = lc.exposure
        grid = np.linspace(-0.3 * P, 0.8 * P, 6000)
        curves = posterior.curves(grid, self.task['epoch'])
        median = np.median(curves, axis=0)
        low, high = np.percentile(curves, [16., 84.], axis=0)
        x, y, err = lc.fold(dt)

        def solve(shift):
            at_data = np.interp(x - shift, grid, median)
            if posterior.nearby:
                coef, chi2, _ = linear_fit(y, err, np.column_stack([at_data, np.ones_like(x)]))
                return float(coef[0]), float(coef[1]), chi2
            coef, chi2, _ = linear_fit(y - at_data, err, np.ones((len(x), 1)))
            return 1., float(coef[0]), chi2

        limit = MAX_SHIFT * self.task['duration']
        shift = float(minimize_scalar(lambda s: solve(s)[2], bounds=(-limit, limit), method='bounded').x)
        scale, baseline, _ = solve(shift)
        self.curve = lambda t: scale * np.interp(t - shift, grid, median) + baseline
        self.band = lambda t: (scale * np.interp(t - shift, grid, low) + baseline, scale * np.interp(t - shift, grid, high) + baseline)
        self.label = 'EB model (posterior)'
        self.t_sec = posterior.secondary_offset() + shift
        self.record.update(route='posterior', scenario=posterior.scenario, scenario_share=posterior.share,
                           dilution_scale=float(scale), shift_min=1440. * shift)


    def from_fit(self):
        """Two-Gaussian fit, then the two-sphere fit started from it (not for contact binaries); the better one is drawn."""
        task, P, lc = self.task, self.period_orb, self.lc
        x, y, err = lc.fold(dt)
        width_1 = task['width_1'] if np.isfinite(task['width_1']) else min(task['duration'] / P, 0.4)
        width_2 = task['width_2'] if np.isfinite(task['width_2']) else width_1
        gauss = TwoGaussianFit(x, y, err, P, task['secondary_phase'] if np.isfinite(task['secondary_phase']) else 0.5, width_1, width_2)
        gauss.run()
        chi2_gauss = float(gauss.minuit.fval)
        self.curve, self.label, self.t_sec = gauss.curve, 'two-Gaussian model', gauss.minuit.values['mu_2'] * P
        self.record.update(route='two-Gaussian fit', twog_chi2=chi2_gauss, twog_depth=float(gauss.coef[1]),
                           twog_depth_err=float(gauss.coef_err[1]), twog_width=gauss.windows()[0][1] / 0.75,
                           twog_shift=gauss.windows()[0][0])
        if self.wide or not self.physical:
            return
        physical = PhysicalFit(lc.star, lc.exposure, x, y, err, P, gauss.windows())
        physical.run(gauss, task['duration'])
        values = tuple(physical.minuit.values)
        chi2 = float(np.sum(((y - physical.curve(x, *values)) / err) ** 2))
        self.record.update(phys_chi2=chi2)
        if chi2 <= PHYSICAL_MARGIN * chi2_gauss:
            self.curve, self.label, self.t_sec = (lambda t: physical.curve(t, *values)), 'EB model (fit)', physical.secondary_time()
            self.record.update(route='physical fit')


    def place_secondary(self):
        """A curve without a visible secondary does not say where to look: the table's secondary phase, else 0.5."""
        task, P = self.task, self.period_orb
        if self.curve is not None:
            _, _, err = self.lc.fold(dt)
            level = np.median(self.curve(np.linspace(-0.25 * P, 0.75 * P, 4000)))
            if level - self.curve(np.array([self.t_sec]))[0] >= 3. * np.median(err):
                self.record['secondary_panel'] = 'model'
                return
        phase_2 = task['secondary_phase'] if np.isfinite(task['secondary_phase']) else 0.5
        self.t_sec = phase_2 * task['period'] if P < 1.5 * task['period'] else task['period']
        self.record['secondary_panel'] = 'table'


    def goodness(self):
        """Per-sample scatter around the curve (scipy's normal-scaled median absolute deviation; 1 would be unit
        variance), and in each eclipse window the depth ratio data / model (least-squares scale of the model's signal)."""
        from scipy.stats import median_abs_deviation
        lc = self.lc
        model = self.curve(lc.phase_time()) if self.curve is not None else np.zeros(len(lc.flux))
        self.noise = float(median_abs_deviation(lc.flux - model, scale='normal'))
        stats = dict(noise_per_sample=self.noise)
        if self.curve is None:
            return stats
        x, y, err = lc.fold(dt)
        at_data = self.curve(x)
        stats['red_chi2'] = float(np.sum(((y - at_data) / (self.noise * err)) ** 2) / max(len(x) - 1, 1))
        for tag, centre in (('primary', 0.), ('secondary', self.t_sec)):
            inside = np.abs(x - centre) < 0.6 * self.task['duration']
            if inside.sum() < 3 or inside.all():
                continue
            signal = at_data[inside] - np.median(at_data[~inside])
            weight = 1. / err[inside] ** 2
            if np.sum(signal ** 2 * weight) > 0:
                stats[tag + '_depth_ratio'] = float(np.sum((y[inside] - np.median(at_data[~inside])) * signal * weight) / np.sum(signal ** 2 * weight))
            stats[tag + '_model_depth'] = float(-signal.min())
        return stats


    def draw(self, png):
        """Detached: primary and secondary zooms on top, the whole orbit below. Contact: the whole orbit only. Data navy
        open circles with the measured scatter as error bars; model black dashed; posterior band skyblue."""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import rc
        from matplotlib.gridspec import GridSpec
        from matplotlib.ticker import MaxNLocator
        rc('mathtext', **{'fontset': 'cm'})
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['DejaVu Serif']
        plt.rcParams['font.size'] = 15

        P, lc = self.period_orb, self.lc
        half = min(max(1.5 * self.task['duration'], 0.15), 0.25 * P)
        coarse, fine = lc.fold(max(P / 200., dt)), lc.fold(dt)
        if self.wide:
            fig = plt.figure(figsize=(9, 5))
            panels = [(fig.add_subplot(1, 1, 1), None)]
        else:
            fig = plt.figure(figsize=(9, 7.6))
            grid = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1.])
            panels = [(fig.add_subplot(grid[0, 0]), 0.), (fig.add_subplot(grid[0, 1]), self.t_sec), (fig.add_subplot(grid[1, :]), None)]

        for ax, centre in panels:
            first = ax is panels[0][0]
            if centre is None:
                x, y, err = coarse
                times, offset, scale = np.linspace(-0.25 * P, 0.75 * P, 3000), 0., P
                ax.set_xlabel('Orbital phase')
                ax.set_xlim(-0.25, 0.75)
            else:
                inside = np.abs(fine[0] - centre) < half
                x, y, err = fine[0][inside], fine[1][inside], fine[2][inside]
                times, offset, scale = centre + np.linspace(-half, half, 1500), centre, 1.
                ax.set_xlabel(r'$t - t_{\rm I}$ [days]' if centre == 0. else r'$t - t_{\rm II}$ [days]')
                ax.set_xlim(-half, half)
                ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.errorbar((x - offset) / scale, y, yerr=self.noise * err, fmt='o', mfc='white', mec='navy', ecolor='navy', ms=5, lw=1.0,
                        capsize=2, label='data' if first else None, zorder=2)
            if self.band is not None:
                low, high = self.band(times)
                ax.fill_between((times - offset) / scale, low, high, color='skyblue', alpha=0.6, lw=0, zorder=1,
                                label='68 % of samples' if first else None)
            if self.curve is not None:
                ax.plot((times - offset) / scale, self.curve(times), color='black', ls='--', lw=2.5, zorder=3,
                        label=self.label if first else None)
            ax.grid(True, alpha=0.25)
        panels[0][0].set_ylabel('Folded Flux (Unit Variance)')
        if not self.wide:
            panels[2][0].set_ylabel('Folded Flux (Unit Variance)')
        handles, labels = panels[0][0].get_legend_handles_labels()
        if self.wide:                                        # the orbit fills the panel: legend above it
            panels[0][0].legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=12, frameon=False)
        else:
            panels[-1][0].legend(handles, labels, loc='lower right', fontsize=12, framealpha=0.9)
        if self.title:
            fig.suptitle(self.title, fontsize=11)
        fig.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(png)), exist_ok=True)
        fig.savefig(png + '.tmp.png', dpi=130)
        plt.close(fig)
        os.replace(png + '.tmp.png', png)


# ----------------------------------------------------------------------------------------------------------------------
# Driver

def png_path(figure):
    return OUT + 'figures/' + figure.replace('plots_eb/', '').replace('/', '_').replace('.webp', '.png')


def one(task):
    """Figure of one row; never raises."""
    warnings.filterwarnings('ignore')
    start = time.time()
    figure = EBFigure(task)
    try:
        record = dict(figure.run(png_path(task['figure'])), error='')
    except Exception as e:
        record = dict(figure.record, error=repr(e)[:300])
    record['seconds'], record['time'] = time.time() - start, time.time()
    return record


def draw(workers, part, only):
    from multiprocessing import get_context
    os.makedirs(OUT + 'figures', exist_ok=True)
    todo = [t for t in EBRows().tasks(only, part) if not os.path.exists(png_path(t['figure']))]
    print(len(todo), 'figures to draw,', workers, 'workers', flush=True)
    start = time.time()
    results = map(one, todo) if workers == 1 else get_context('spawn').Pool(workers, maxtasksperchild=20).imap_unordered(one, todo, chunksize=1)
    records = OUT + 'records_%d_of_%d.csv' % part
    for i, record in enumerate(results):
        pd.DataFrame([record]).to_json(records, orient='records', lines=True, mode='a')
        if record['error'] or (i + 1) % 200 == 0 or workers == 1:
            print('[%d/%d] %s %s %.0f s %s (%.1f min)' % (i + 1, len(todo), record['name'], record.get('route', ''), record['seconds'],
                                                         record['error'], (time.time() - start) / 60), flush=True)
    print('done in %.1f min' % ((time.time() - start) / 60), flush=True)


def collect():
    """The last record of every row, with whether its figure file exists."""
    import glob
    frames = [pd.read_json(path, orient='records', lines=True) for path in glob.glob(OUT + 'records_*.csv')]
    records = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=['figure'])
    if len(records):
        if 'time' not in records:
            records['time'] = np.nan
        records = records.sort_values('time', na_position='first', kind='stable').drop_duplicates('figure', keep='last')
        records['drawn'] = [os.path.exists(png_path(f)) for f in records['figure']]
    return records


def summary():
    rows, records = EBRows().table, collect()
    print('rows %d, on archive stars %d; records %d, drawn %d' % (len(rows), (rows['in_archive'] == True).sum(), len(records),
                                                                 records['drawn'].sum() if len(records) else 0))
    if not len(records):
        return
    done = records[records['drawn']]
    print(done.groupby('route').agg(figures=('figure', 'size'), seconds=('seconds', 'median'), noise=('noise_per_sample', 'median'),
                                    red_chi2=('red_chi2', 'median'), primary_depth_ratio=('primary_depth_ratio', 'median')).round(3).to_string())
    posterior = done[done['route'] == 'posterior']
    if len(posterior):
        print('posterior curves, primary depth data / model by scenario (5, 50, 95 %):')
        print(posterior.groupby('scenario')['primary_depth_ratio'].quantile([.05, .5, .95]).unstack().round(2).to_string())
    errors = records['error'].fillna('').astype(str)
    print('errors:', errors[errors != ''].str[:80].value_counts().head(8).to_dict())


def install(apply):
    records = collect()
    done = records[records['drawn']] if len(records) else records
    print('figures drawn %d of %d records; errors %d' % (len(done), len(records), (records['error'].fillna('').astype(str) != '').sum() if len(records) else 0))
    if not apply:
        print('dry run; pass --apply to write the WebP files and eb_figures.tsv')
        return
    written = 0
    for r in done.itertuples():                              # only figures drawn since their WebP was written
        source, target = png_path(r.figure), os.path.join(HERE, r.figure)
        if not os.path.exists(target) or os.path.getmtime(target) < os.path.getmtime(source):
            site_figures.to_webp(source, target)
            written += 1
    print('WebP files written:', written)
    done.drop(columns='drawn').to_csv(RESULTS, sep='\t', index=False)
    print('installed', len(done), 'figures; now run build_verified_catalogs.py')


if __name__ == '__main__':
    stage = sys.argv[1] if len(sys.argv) > 1 else ''
    if stage == 'list':
        table = EBRows().table
        print(table.groupby(['route', table['in_archive'].map({True: 'in the archive', False: 'not in the archive'})]).size().unstack(fill_value=0).to_string())
    elif stage == 'draw':
        args = sys.argv[3:]
        part = tuple(int(v) for v in args.pop(0).split('/')) if args and '/' in args[0] else (0, 1)
        draw(int(sys.argv[2]) if len(sys.argv) > 2 else 4, part, [int(a) for a in args] or None)
    elif stage == 'summary':
        summary()
    elif stage == 'install':
        install('--apply' in sys.argv)
    else:
        sys.exit(__doc__)
