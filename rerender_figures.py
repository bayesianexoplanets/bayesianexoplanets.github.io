"""rerender_figures.py -- redraw every transit figure of the website with larger type, folded at the listed ephemeris.

`tess_planets.csv` rows are folded at their listed (verified) period and epoch with the fit's half-duration, the new
candidates (`tois_new.csv`) at their own Period / Phase / Tau, both by `pipeline.pretty_plot.pretty_transit_plot` on the
website's SPOC light curves. Only the type size differs from the pipeline's figures (FONT_SIZE, the house style).

    python rerender_figures.py render [workers] [TIC ...]   PNGs to scratch, resumable (an existing PNG is skipped)
    python rerender_figures.py install [--apply]            PNG -> JPG into plots/ and plots_new/ (dry run by default)
"""
import os
import sys
import time
import warnings

for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_var, '1')
os.environ.pop('TESS_SHARD_ARCHIVE', None)             # the website's light curves

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
from constants import scratch

OUT = scratch + 'Rerender_20261120/'                     # Rerender_20261119/: drawn with the SPOC files' grid read as 600 s
SUMMARY = OUT + 'summary.csv'
FONT_SIZE = 15                                         # the house style; the pipeline's figures use 11
LEGACY_STEP = 0.001388882 * 5                          # days: the grid the website's SPOC files were built on


def tasks(only=None):
    """(jpg path relative to the site, TIC, period, epoch or phase, epoch is BTJD?, tau, SNR) of every figure."""
    planets = pd.read_csv(os.path.join(HERE, 'tess_planets.csv'), low_memory=False)
    planets = planets[planets['Data'] == 'SPOC']           # archive-fitted rows are drawn by known_figures.py
    new = pd.read_csv(os.path.join(HERE, 'tois_new.csv'), low_memory=False)
    rows = [(r.plot, int(r.TIC), float(r.Period), float(r.Epoch), True, float(r.Tau_fit), float(r.SNR)) for r in planets.itertuples()]
    rows += [('plots_new/%d/%d.jpg' % (r.TIC, r.cand_idx), int(r.TIC), float(r.Period), float(r.Phase), False, float(r.Tau), float(r.SNR))
             for r in new.rename(columns={'_cand_idx': 'cand_idx'}).itertuples()]
    rows = [t for t in rows if np.all(np.isfinite(t[2:4])) and np.isfinite(t[5]) and t[5] > 0]
    return [t for t in rows if only is None or t[1] in only]


def png_path(plot):
    return OUT + os.path.splitext(plot)[0].replace('/', '_') + '.png'


def grid_stretch():
    """Grid time per day of BTJD on the light curves the loaders currently read.

    The pipeline counts time as t_start + Indices * dt. The all-sky archive was built with that step; the website's
    SPOC files were built with LEGACY_STEP, so on them a time elapsed since t_start reads longer by dt / LEGACY_STEP.
    """
    from constants import dt
    return 1. if os.environ.get('TESS_SHARD_ARCHIVE') else dt / LEGACY_STEP


def draw(png, tic, period, when, is_epoch, tau, snr):
    """Render one figure to `png` on the light curves the loaders currently read; never raises."""
    warnings.filterwarnings('ignore')
    import matplotlib
    matplotlib.use('Agg')
    t0 = time.time()
    try:
        from pipeline.pretty_plot import pretty_transit_plot
        from TESS.load_tess import load_stellar_data
        import matplotlib.pyplot as plt
        plt.rcParams['font.size'] = FONT_SIZE
        stretch = grid_stretch()
        elapsed = when - float(load_stellar_data(tic)[0]) if is_epoch else when
        phase = float(np.mod(elapsed * stretch, period * stretch))
        summary = pretty_transit_plot(tic, period * stretch, phase, tau * stretch, png, snr=snr)
        return dict(TIC=tic, model=str(summary.get('model', 'none')), seconds=time.time() - t0, error='')
    except Exception as e:
        return dict(TIC=tic, model='none', seconds=time.time() - t0, error=repr(e)[:200])


def one(task):
    plot, tic, period, when, is_epoch, tau, snr = task
    return dict(plot=plot, **draw(png_path(plot), tic, period, when, is_epoch, tau, snr))


def render(workers, only=None):
    from multiprocessing import get_context
    os.makedirs(OUT, exist_ok=True)
    todo = [t for t in tasks(only) if not os.path.exists(png_path(t[0]))]
    print(len(todo), 'figures to render,', workers, 'workers', flush=True)
    t_start = time.time()
    results = map(one, todo) if workers == 1 else get_context('spawn').Pool(workers, maxtasksperchild=20).imap_unordered(one, todo, chunksize=1)
    for i, result in enumerate(results):
        pd.DataFrame([result]).to_csv(SUMMARY, mode='a', header=not os.path.exists(SUMMARY), index=False)
        if (i + 1) % 200 == 0 or result['error'] or workers == 1:
            print('[%d/%d] %s %.0f s %s (%.1f min)' % (i + 1, len(todo), result['plot'], result['seconds'], result['error'], (time.time() - t_start) / 60), flush=True)
    print('done in %.1f min' % ((time.time() - t_start) / 60), flush=True)


def install(apply):
    import site_figures
    every = tasks()
    ready = [t for t in every if os.path.exists(png_path(t[0]))]
    print('%d of %d figures rendered (planets %d, new candidates %d)'
          % (len(ready), len(every), sum(t[0].startswith('plots/') for t in ready), sum(t[0].startswith('plots_new/') for t in ready)))
    if not apply:
        print('dry run; pass --apply to write the JPGs')
        return
    for t in ready:
        site_figures.install(png_path(t[0]), HERE, t[0])
    print('installed', len(ready), '; bump PLOT_VERSION in index.html')


if __name__ == '__main__':
    stage = sys.argv[1]
    if stage == 'render':
        only = {int(a) for a in sys.argv[3:]} or None
        render(int(sys.argv[2]) if len(sys.argv) > 2 else 4, only)
    elif stage == 'install':
        install('--apply' in sys.argv)
    else:
        sys.exit(__doc__)
