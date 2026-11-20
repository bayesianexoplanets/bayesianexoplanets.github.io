"""true_time_catalog.py -- the planet fits on the samples' true times, as the website lists them.

Since 2026-11-20 the known-object fit runs on the true BTJD of every sample instead of the 10-min grid (logs/DESIGN.md).
Every planet row of the site was refitted that way (allsky_search/true_time_fits_20261120/: three iterations, then a
check against independently measured BTJD mid-times, `tables/decision_full.tsv`).

    python true_time_catalog.py table                  true_time_fits.tsv: every refitted row and whether the site adopts it
    python true_time_catalog.py render WORKERS [TIC ...]   figures of the rows on their true times, to scratch (resumable)
    python true_time_catalog.py install [--apply]      the true-time figures of the adopted rows -> plots/, plots_known/

A row is adopted when the mid-times agree with the fit; or when it has no usable mid-times but the fit converged, lost
no SNR and stayed on the listed signal; or when its mid-times only fail to follow a line while the fit is healthy (the
linear ephemeris is listed until the separate TTV analysis); and its true-time figure exists. For an adopted row
`build_verified_catalogs.py` lists the fit's period, epoch, errors, SNR, radius and transit count, and
`site_figures.install` lets only a true-time fold replace its figure.
"""
import glob
import os
import sys
import time
import warnings

for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from constants import scratch
import site_figures

RERUN = scratch + 'allsky_search/true_time_fits_20261120/'
OUT = scratch + 'Rerender_truetime_20261120/'
ARCHIVE = scratch + 'allsky/archive'
TABLE = os.path.join(HERE, site_figures.TRUE_TIME)
FONT_SIZE = 15                                          # the house style, as rerender_figures
FLAGGED = 'ttv: redo separately'


def png_path(plot):
    return OUT + os.path.splitext(plot)[0].replace('/', '_') + '.png'


def last_iterations():
    """The last iteration of every refitted row (its phase, errors, SNR, radius and transit count)."""
    frames = [pd.read_csv(path, sep='\t').assign(data=data) for data in ('spoc', 'archive')
              for path in sorted(glob.glob(RERUN + 'full/%s/*.tsv' % data)) if not path.endswith('_check.tsv')]
    every = pd.concat(frames, ignore_index=True)
    return every[every['iteration'] > 0].sort_values('iteration').groupby(['data', 'key']).last().reset_index()


def table():
    decision = pd.read_csv(RERUN + 'tables/decision_full.tsv', sep='\t', low_memory=False)
    fit = last_iterations()[['data', 'key', 'phase', 't_start', 'SNR', 'err_tau', 'radius', 'n_transits']].rename(columns={'n_transits': 'n_fit'})
    rows = decision.merge(fit, on=['data', 'key'], how='left', validate='one_to_one')
    assert np.allclose(rows['n_fit'], rows['n_transits'], equal_nan=True), 'the two tables disagree on a row'
    healthy = rows['converged'] & rows['same_signal'] & (rows['snr_mf_true'] >= rows['snr_mf_listed'] - 0.1)   # as decision_full.py
    wanted = rows['action'].str.startswith('adopt') | ((rows['action'] == FLAGGED) & healthy)
    finite = np.isfinite(rows[['P_true', 'E_true', 'tau_true', 'SNR', 'phase']]).all(axis=1) & (rows['tau_true'] > 0)
    drawn = np.array([os.path.exists(png_path(key)) for key in rows['key']])
    out = pd.DataFrame({
        'plot': rows['key'], 'tic': rows['TIC'], 'data': rows['data'], 'action': rows['action'],
        'adopted': wanted & finite & drawn, 'nonlinear': rows['nonlinear'], 'gate': rows['gate'].fillna(''),
        'period': rows['P_true'], 'epoch': rows['E_true'], 'tau': rows['tau_true'], 'phase': rows['phase'], 'SNR': rows['SNR'],
        'err_period': rows['err_P'], 'err_epoch': rows['err_E'], 'err_tau': rows['err_tau'], 'radius': rows['radius'],
        'n_transits': rows['n_transits'], 'period_before': rows['P_listed'], 'epoch_before': rows['E_listed'],
        'tau_before': rows['tau_listed'], 'snr_before': rows['snr_mf_listed'], 'snr_after': rows['snr_mf_true'],
        'chi2_line': rows['chi2_line'], 'chi2_before': rows['chi2_listed'], 'chi2_after': rows['chi2_true'], 'dof': rows['dof']})
    out.to_csv(TABLE + '.tmp', sep='\t', index=False)
    os.replace(TABLE + '.tmp', TABLE)
    print('%s: %d rows, adopted %d (wanted %d, of these without a figure %d, not finite %d)'
          % (os.path.basename(TABLE), len(out), out['adopted'].sum(), wanted.sum(), (wanted & ~drawn).sum(), (wanted & ~finite).sum()))
    print(pd.crosstab(out['action'], out['adopted']).to_string())


def draw_one(task):
    """The figure of one row at its true-time ephemeris, folded on the samples' true times; never raises."""
    data, tic, plot, period, epoch, tau, snr = task
    warnings.filterwarnings('ignore')
    import matplotlib
    matplotlib.use('Agg')
    began = time.time()
    out = dict(data=data, TIC=tic, key=plot, error='')
    try:
        if data == 'archive':
            os.environ['TESS_SHARD_ARCHIVE'] = ARCHIVE
        else:
            os.environ.pop('TESS_SHARD_ARCHIVE', None)             # the website's SPOC files
        from pipeline.pretty_plot import pretty_transit_plot
        from TESS.load_tess import load_stellar_data
        import matplotlib.pyplot as plt
        plt.rcParams['font.size'] = FONT_SIZE
        t_start = float(load_stellar_data(tic)[0])
        summary = pretty_transit_plot(tic, period, float(np.mod(epoch - t_start, period)), tau, png_path(plot), snr=snr, true_time=True)
        out.update({k: summary.get(k, np.nan) for k in ('fold_snr', 'T14', 'rp', 't0', 'n_epochs', 'true_time', 'model')})
    except Exception as e:
        out['error'] = repr(e)[:200]
    out['seconds'] = time.time() - began
    return out


def render(workers, only=None):
    """Draw the figure of every row the table would adopt (it is adopted only once its figure exists)."""
    from multiprocessing import get_context
    os.makedirs(OUT, exist_ok=True)
    decision = pd.read_csv(RERUN + 'tables/decision_full.tsv', sep='\t', low_memory=False)
    wanted = decision[~decision['action'].str.startswith('weak') & np.isfinite(decision['P_true']) & (decision['tau_true'] > 0)]
    todo = [(r.data, int(r.TIC), r.key, float(r.P_true), float(r.E_true), float(r.tau_true), float(r.snr_mf_true))
            for r in wanted.itertuples() if not os.path.exists(png_path(r.key)) and (only is None or int(r.TIC) in only)]
    print(len(todo), 'figures to draw of', len(wanted), 'rows,', workers, 'workers', flush=True)
    began = time.time()
    results = map(draw_one, todo) if workers == 1 else get_context('spawn').Pool(workers, maxtasksperchild=10).imap_unordered(draw_one, todo, chunksize=1)
    for i, result in enumerate(results):
        log = OUT + 'summary_site.tsv'
        pd.DataFrame([result]).to_csv(log, sep='\t', mode='a', header=not os.path.exists(log), index=False)
        if (i + 1) % 200 == 0 or result['error'] or workers == 1:
            print('[%d/%d] %s %.0f s %s (%.1f min)' % (i + 1, len(todo), result['key'], result['seconds'], result['error'], (time.time() - began) / 60), flush=True)


def install(apply):
    rows = pd.read_csv(TABLE, sep='\t')
    adopted = rows[rows['adopted']]
    ready = [plot for plot in adopted['plot'] if os.path.exists(png_path(plot))]
    print('%d adopted rows, %d with a figure (SPOC %d, archive %d)' % (len(adopted), len(ready), sum(p.startswith('plots/') for p in ready),
                                                                        sum(p.startswith('plots_known/') for p in ready)))
    if not apply:
        print('dry run; pass --apply to write the figures')
        return
    for plot in ready:
        site_figures.install(png_path(plot), HERE, plot, true_time=True)
    print('installed', len(ready), '; bump PLOT_VERSION in index.html and run build_verified_catalogs.py')


if __name__ == '__main__':
    stage = sys.argv[1] if len(sys.argv) > 1 else ''
    if stage == 'table':
        table()
    elif stage == 'render':
        render(int(sys.argv[2]) if len(sys.argv) > 2 else 4, {int(a) for a in sys.argv[3:]} or None)
    elif stage == 'install':
        install('--apply' in sys.argv)
    else:
        sys.exit(__doc__)
