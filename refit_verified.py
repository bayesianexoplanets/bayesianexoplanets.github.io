"""refit_verified.py -- refit the planet rows whose fit is stale at their verified ephemeris, and redraw their figures.

`tess_planets.csv` lists the verified period, epoch and full duration; for `fit_stale` rows the SNR, radius, p-value and
figure still belong to the earlier fit. Here each such star goes through the known-object stage of the pipeline on the
website's own light curves (the downloaded SPOC files, the data its null searches ran on), every catalogued signal of
the star seeded at its verified ephemeris with the smallest ephemeris box. The figure is then folded at the LISTED
ephemeris with the fitted half-duration, so table and figure agree. The null distributions are not redone.

    python refit_verified.py list                 the rows to refit
    python refit_verified.py fit [workers]        fits and figures to scratch, resumable (a finished star is skipped)
    python refit_verified.py install [--apply]    figures into plots/, results into refits.tsv (dry run by default);
                                                  then run build_verified_catalogs.py
"""
import os
import sys
import time
import warnings

for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(_var, '1')
os.environ.pop('TESS_SHARD_ARCHIVE', None)             # the website's light curves, not the all-sky archive

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, REPO)
from constants import scratch

FOLDER = 'Refit_20261119'                              # under scratch
OUT = scratch + FOLDER + '/'
HYPERPARAMS = dict(do_nonstat_mf=False, test_recovery=False, cutoff_fa_eliminate=8., period_density=300000.0 / 2.0,
                   show=False, snrd_by_origin=True, cutoff_num=10)       # submissions/run_pipeline.py, batch 0, no plots
TIGHT_PERIOD_ERR = 1e-7                                # relative: the ephemeris box keeps its minimum size
EARTH_PER_SUN = 0.009167888                            # Earth radius in solar radii
AGREE = 0.5                                            # listed durations between fit and listed transit times


def catalogs():
    planets = pd.read_csv(os.path.join(HERE, 'tess_planets.csv'), low_memory=False)
    planets = planets[planets['Data'] == 'SPOC']           # the rows this pipeline fitted on the SPOC files
    ebs = pd.read_csv(os.path.join(HERE, 'tess_ebs.csv'), low_memory=False)
    return planets, ebs


def star_table(tic, planets, ebs):
    """The star's catalogued signals as a known-object table of the search (`pipeline.search_masks.KnownObjects`)."""
    rows = []
    for r in planets[planets['TIC'] == tic].itertuples():
        ror = EARTH_PER_SUN * r.Radius_planet / r.Radius if np.isfinite(r.Radius_planet) and np.isfinite(r.Radius) and r.Radius > 0 else np.nan
        rows.append(dict(name=r.plot, kind='planet', object_class='planet candidate', eb_flag=False, period=r.Period, epoch=r.Epoch,
                         duration=r.Duration, ror=ror, secondary_phase=np.nan, secondary_duration=np.nan, toi=r.TOI))
    for r in ebs[ebs['TIC'] == tic].itertuples():
        rows.append(dict(name='EB ' + str(r.Name), kind='eb', object_class='EB', eb_flag=True, period=r.Period, epoch=r.Epoch,
                         duration=r.Duration, ror=np.nan, secondary_phase=r.Secondary_phase, secondary_duration=r.Secondary_duration, toi=np.nan))
    table = pd.DataFrame(rows)
    table = table[np.isfinite(table['period']) & (table['period'] > 0) & np.isfinite(table['epoch'])]
    return table.assign(tic=tic, masked=True, frozen=False, source='website', depth_ppm=np.nan, ttv_widen=np.nan,
                        period_err=TIGHT_PERIOD_ERR * table['period'], epoch_err=0.)


def transit_gap(period, epoch, fit_period, fit_epoch, times):
    """Largest distance, days, between the listed transits that fall on observed samples and the fitted ones."""
    numbers = np.unique(np.round((times - epoch) / period))
    centres = epoch + numbers * period
    centres = centres[np.min(np.abs(times[None, :] - centres[:, None]), axis=1) < 0.02]   # a sample within half an hour
    if len(centres) == 0:
        return np.nan
    return float(np.max(np.abs(np.mod(centres - fit_epoch + 0.5 * fit_period, fit_period) - 0.5 * fit_period)))


def fit_record(tic, info, listed, t_start, observed):
    """One fitted object as a table row: its fit, its listed ephemeris and the distance between their transit times."""
    period, phase, tau = (float(x) for x in info.params)
    cov = np.asarray(getattr(info, 'cov', np.full(6, np.nan)), dtype=float)
    radius = np.asarray(getattr(info, 'radius', np.full(3, np.nan)), dtype=float)
    return dict(TIC=tic, listed_period=listed['period'], listed_epoch=listed['epoch'], listed_duration=listed['duration'],
                period=period, epoch=t_start + phase, phase=phase, tau=tau, SNR=float(info.snr),
                err_period=np.sqrt(cov[0]) if cov[0] > 0 else np.nan, err_epoch=np.sqrt(cov[1]) if cov[1] > 0 else np.nan,
                err_tau=np.sqrt(cov[2]) if cov[2] > 0 else np.nan, amplitude=float(getattr(info, 'amplitude', np.nan)),
                radius=radius[0], radius_errp=radius[1], radius_errm=radius[2],
                n_transits=int(getattr(info, 'num_valid_transits', -1)), has_ttvs=hasattr(info, 'ttvs'),
                gap=transit_gap(listed['period'], listed['epoch'], period, t_start + phase, observed))


def one(tic):
    """Fit one star and draw the figures of its stale rows; never raises."""
    warnings.filterwarnings('ignore')
    import matplotlib
    matplotlib.use('Agg')
    t0 = time.time()
    try:
        from constants import dt
        from pipeline.main import mainn
        from pipeline.hyperparameters import Hyperparams
        from pipeline.search_masks import SearchMasks, KnownObjects
        from pipeline.pretty_plot import pretty_transit_plot
        from TESS.load_tess import load_blocks
        planets, ebs = catalogs()
        table = star_table(tic, planets, ebs).sort_values('tic').reset_index(drop=True)
        path = OUT + 'tables/%d.tsv' % tic
        table.to_csv(path, sep='\t', index=False)
        masks = SearchMasks(KnownObjects(path), fit_known=True)
        infos = mainn(tic, FOLDER, Hyperparams(**HYPERPARAMS), use_planet_data=False, masks=masks, known_planets_only=True)

        Flux, invVar, Indices, upsampled, stellar = load_blocks(tic)
        t_start = float(stellar[0])
        observed = (t_start + Indices * dt)[invVar > 0.5]
        stale = set(planets.loc[(planets['TIC'] == tic) & planets['fit_stale'], 'plot'])
        out = []
        for info in infos:
            row = table.iloc[info.catalog['position']]
            result = dict(plot=row['name'], **fit_record(tic, info, row, t_start, observed), figure='', error='')
            if row['name'] in stale and np.isfinite(result['SNR']):
                png = OUT + 'figures/' + row['name'].replace('/', '_').replace('.jpg', '.png')
                try:
                    pretty_transit_plot(tic, float(row['period']), float(np.mod(row['epoch'] - t_start, row['period'])), result['tau'], png, snr=result['SNR'])
                    result['figure'] = png
                except Exception as e:
                    result['error'] = 'figure: ' + repr(e)[:150]
            out.append(result)
        missing = set(table['name']) - {r['plot'] for r in out}
        out += [dict(TIC=tic, plot=name, error='not fitted') for name in missing]
        frame = pd.DataFrame(out)
    except Exception as e:
        frame = pd.DataFrame([dict(TIC=tic, plot='', error=repr(e)[:300])])
    frame['seconds'] = time.time() - t0
    frame.to_csv(OUT + 'fits/%d.csv.tmp' % tic, index=False)
    os.replace(OUT + 'fits/%d.csv.tmp' % tic, OUT + 'fits/%d.csv' % tic)
    return tic, len(frame), float(time.time() - t0), '; '.join(sorted(set(frame['error'].dropna().astype(str)) - {''}))


def stars():
    planets, _ = catalogs()
    return sorted(int(t) for t in planets.loc[planets['fit_stale'], 'TIC'].unique())


def fit(workers, only=None):
    from multiprocessing import get_context
    for sub in ('tables', 'fits', 'figures'):
        os.makedirs(OUT + sub, exist_ok=True)
    todo = [t for t in (only or stars()) if not os.path.exists(OUT + 'fits/%d.csv' % t)]
    print(len(todo), 'stars to fit,', workers, 'workers', flush=True)
    t_start = time.time()
    if workers == 1:
        results = map(one, todo)
    else:
        pool = get_context('spawn').Pool(workers, maxtasksperchild=10)
        results = pool.imap_unordered(one, todo, chunksize=1)
    for i, (tic, n, seconds, error) in enumerate(results):
        if error or (i + 1) % 25 == 0 or workers == 1:
            print('[%d/%d] TIC %d: %d rows, %.0f s %s (%.1f min)' % (i + 1, len(todo), tic, n, seconds, error, (time.time() - t_start) / 60), flush=True)
    print('done in %.1f min' % ((time.time() - t_start) / 60), flush=True)


def collect():
    """Every fitted row of the stale stars, with whether the fit stayed at the listed ephemeris."""
    import glob
    frames = [pd.read_csv(p) for p in glob.glob(OUT + 'fits/*.csv')]
    fits = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    if len(fits):
        fits['agrees'] = fits['gap'] <= AGREE * fits['listed_duration']
    return fits


def install(apply):
    import site_figures
    site = pd.read_csv(os.path.join(HERE, 'tois.csv'), low_memory=False, usecols=['TIC', 'SNR'])
    site['plot'] = 'plots/' + site['TIC'].astype(str) + '/' + site.groupby('TIC').cumcount().astype(str) + '.jpg'
    fits = collect()
    targets = fits[fits['figure'].fillna('') != '']                        # a figure was drawn for every row that was stale
    merged = targets.merge(site.rename(columns={'SNR': 'SNR_before'})[['plot', 'SNR_before']], on='plot', how='left')
    print('refitted rows %d: fit agrees with the listed ephemeris %d, errors %d'
          % (len(merged), (merged['agrees'] == True).sum(), (merged['error'].fillna('') != '').sum()))
    print('SNR after / before, rows that agree: quantiles',
          (merged['SNR'] / merged['SNR_before'])[merged['agrees'] == True].quantile([.05, .5, .95]).round(3).to_dict())
    if not apply:
        print('dry run; pass --apply to write refits.tsv and the figures')
        return
    merged.to_csv(os.path.join(HERE, 'refits.tsv'), sep='\t', index=False)
    for r in merged.itertuples():
        site_figures.install(r.figure, HERE, r.plot)
    print('installed', len(merged), 'rows; now run build_verified_catalogs.py')


if __name__ == '__main__':
    stage = sys.argv[1]
    if stage == 'list':
        planets, _ = catalogs()
        print(planets.loc[planets['fit_stale'], ['TIC', 'TOI', 'Period', 'Period_fit', 'SNR', 'Ephemeris']].to_string(index=False))
        print(len(stars()), 'stars')
    elif stage == 'fit':
        only = [int(a) for a in sys.argv[3:]] or None
        fit(int(sys.argv[2]) if len(sys.argv) > 2 else 4, only)
    elif stage == 'install':
        install('--apply' in sys.argv)
    else:
        sys.exit(__doc__)
