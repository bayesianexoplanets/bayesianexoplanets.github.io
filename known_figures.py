"""known_figures.py -- fits and figures of the verified objects that have no figure on the website, on the all-sky archive.

The website's figures cover the objects its own pipeline fitted on the downloaded SPOC light curves. Every other
verified planet, planet candidate and eclipsing binary that the all-sky search fits before masking it (a table row that
is masked and not frozen, on a star of the archive) goes through that same known-object stage here: the search's table
and hyperparameters, on the archive's light curves. Its figure is folded at the LISTED ephemeris with the fitted
half-duration, as on the rest of the website.

    python known_figures.py list                          the objects without a figure, by class
    python known_figures.py fit [workers] [TIC ...]       known-object stage per star, resumable (a finished star is skipped)
    python known_figures.py draw [workers] [class ...]    figures of the fitted objects, resumable; classes: planet,
                                                          candidate, eb (default: planet candidate)
    python known_figures.py summary                       fits and figures so far
    python known_figures.py install [--apply]             planets and candidates into the website: `known_fits.tsv` and
                                                          the figures under plots_known/ (dry run by default); then run
                                                          build_verified_catalogs.py

With `--nea` anywhere on the command line every stage works on a second group instead: the confirmed planets without a
TOI number (NASA Exoplanet Archive members that the search's table keeps `frozen`, i.e. masks at their catalog values
without fitting). For their figures they are fitted here like the others, through a copy of the table in which they
are not frozen; results go to `known_fits_nea.tsv`.

Runs in the pipeline's environment (`cosmology`, Python 3.8).
"""
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
import refit_verified as rv
import rerender_figures as rr

ARCHIVE = scratch + 'allsky/archive'
TABLE = scratch + 'allsky_search/known_objects/known_objects_20261114_verified_v2.tsv'
FOLDER = 'Known_figures_20261119'                      # under scratch
OUT = scratch + FOLDER + '/'
TARGETS = OUT + 'targets.tsv'
CLASSES = {'planet': 'planet', 'candidate': 'planet candidate', 'eb': 'EB'}
KNOWN_FITS = os.path.join(HERE, 'known_fits.tsv')      # read by build_verified_catalogs.py
_state = {}


def use_archive():
    """Point the loaders at the all-sky archive (the website scripts imported above point them at the SPOC files)."""
    os.environ['TESS_SHARD_ARCHIVE'] = ARCHIVE


def nea():
    """Whether this run works on the confirmed planets without a TOI number (`--nea`); inherited by the worker processes."""
    return os.environ.get('KNOWN_FIGURES_GROUP', '') == 'nea'


def place(kind):
    """Where the current group keeps its fits, target list, table, figure log and website table."""
    tag = '_nea' if nea() else ''
    return {'fits': OUT + 'fits%s/' % tag, 'targets': OUT + 'targets%s.tsv' % tag, 'log': OUT + 'figures%s.csv' % tag,
            'table': OUT + 'table_nea_unfrozen.tsv' if nea() else TABLE, 'site': os.path.join(HERE, 'known_fits%s.tsv' % tag)}[kind]


def targets():
    """The verified objects without a website figure, one row each, with the name of their figure: those the search fits,
    or with `--nea` the frozen planets on archive stars that no website row owns (`build_verified_catalogs.nea_rows`)."""
    from build_verified_catalogs import VerifiedCatalogs
    from analysis import known_objects as ko
    catalogs = VerifiedCatalogs()
    table, site = catalogs.table, catalogs.site[catalogs.site['row'] >= 0]
    assert not table.duplicated(['tic', 'name']).any(), 'a star lists one name twice'
    own = np.abs(site['Period'].to_numpy() / table.loc[site['row'], 'period'].to_numpy() - 1.) < ko.SAME_PERIOD
    shown = set(site.loc[own, 'row'])

    usable = table['masked'] & table['in_archive'] & np.isfinite(table['period']) & (table['period'] > 0) & np.isfinite(table['epoch'])
    if nea():
        wanted = table[usable & table['frozen'] & table['object_class'].isin([CLASSES['planet'], CLASSES['candidate']])
                       & ~table.index.isin(set(site['row']))].copy()
        if not os.path.exists(place('table')):                 # the search's table with these planets open to the fit
            os.makedirs(OUT, exist_ok=True)
            pd.read_csv(TABLE, sep='\t', low_memory=False).assign(frozen=lambda t: t['frozen'] & ~t.index.isin(wanted.index)).to_csv(
                place('table'), sep='\t', index=False)
    else:
        wanted = table[usable & ~table['frozen'] & table['object_class'].isin(CLASSES.values()) & ~table.index.isin(shown)].copy()
    wanted['plot'] = 'plots_known/' + wanted['tic'].astype(str) + '/' + table.groupby('tic').cumcount()[wanted.index].astype(str) + '.jpg'
    return wanted[['tic', 'name', 'object_class', 'period', 'epoch', 'duration', 'secondary_phase', 'secondary_duration',
                   'review_signal', 'toi', 'plot']]


def search_masks():
    """The search's own masks, loaded once per process."""
    if 'masks' not in _state:
        from pipeline.search_masks import SearchMasks, KnownObjects
        _state['masks'] = SearchMasks(KnownObjects(place('table')), fit_known=True)
    return _state['masks']


def fit_one(tic):
    """Known-object stage of one star; never raises."""
    warnings.filterwarnings('ignore')
    import matplotlib
    matplotlib.use('Agg')
    use_archive()
    t0 = time.time()
    try:
        from constants import dt
        from pipeline.main import mainn
        from pipeline.hyperparameters import Hyperparams
        from TESS.load_tess import load_blocks
        masks = search_masks()
        infos = mainn(tic, FOLDER, Hyperparams(**rv.HYPERPARAMS), use_planet_data=False, masks=masks, known_planets_only=True)

        Flux, invVar, Indices, upsampled, stellar = load_blocks(tic)
        t_start = float(stellar[0])
        observed = (t_start + Indices * dt)[invVar > 0.5]
        objects = masks.known.objects(tic)
        out = [dict(name=objects.iloc[info.catalog['position']]['name'],
                    **rv.fit_record(tic, info, objects.iloc[info.catalog['position']], t_start, observed),
                    star_radius=float(stellar[7]), error='') for info in infos]
        out += [dict(TIC=tic, name=name, error='not fitted') for name in set(objects['name']) - {r['name'] for r in out}]
        frame = pd.DataFrame(out)
    except Exception as e:
        frame = pd.DataFrame([dict(TIC=tic, name='', error=repr(e)[:300])])
    frame['seconds'] = time.time() - t0
    frame.to_csv(place('fits') + '%d.csv.tmp' % tic, index=False)
    os.replace(place('fits') + '%d.csv.tmp' % tic, place('fits') + '%d.csv' % tic)
    return tic, len(frame), float(time.time() - t0), '; '.join(sorted(set(frame['error'].dropna().astype(str)) - {''}))


def draw_one(task):
    """Figure of one fitted object at its listed ephemeris; never raises."""
    use_archive()
    plot, tic, period, epoch, tau, snr = task
    return dict(plot=plot, **rr.draw(png_path(plot), tic, period, epoch, True, tau, snr))


def png_path(plot):
    return OUT + 'figures/' + plot.replace('plots_known/', '').replace('/', '_').replace('.jpg', '.png')


def pool_map(worker, todo, workers):
    from multiprocessing import get_context
    return map(worker, todo) if workers == 1 else get_context('spawn').Pool(workers, maxtasksperchild=10).imap_unordered(worker, todo, chunksize=1)


def fit(workers, only=None):
    for folder in (place('fits'), OUT + 'figures'):
        os.makedirs(folder, exist_ok=True)
    wanted = targets()
    wanted.to_csv(place('targets'), sep='\t', index=False)
    stars = only or sorted(int(t) for t in wanted['tic'].unique())
    todo = [t for t in stars if not os.path.exists(place('fits') + '%d.csv' % t)]
    print(len(todo), 'stars to fit,', workers, 'workers', flush=True)
    t_start = time.time()
    for i, (tic, n, seconds, error) in enumerate(pool_map(fit_one, todo, workers)):
        if error or (i + 1) % 100 == 0 or workers == 1:
            print('[%d/%d] TIC %d: %d rows, %.0f s %s (%.1f min)' % (i + 1, len(todo), tic, n, seconds, error, (time.time() - t_start) / 60), flush=True)
    print('done in %.1f min' % ((time.time() - t_start) / 60), flush=True)


def collect():
    """The targets with their fit, where one exists, and whether it stayed at the listed ephemeris."""
    import glob
    wanted = pd.read_csv(place('targets'), sep='\t')
    frames = [pd.read_csv(p) for p in glob.glob(place('fits') + '*.csv')]
    fits = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=['TIC', 'name'])
    fits['name'] = fits['name'].fillna('').astype(str)
    merged = wanted.merge(fits.rename(columns={'TIC': 'tic', 'period': 'fit_period', 'epoch': 'fit_epoch'}), on=['tic', 'name'], how='left')
    if 'gap' in merged:
        merged['agrees'] = merged['gap'] <= rv.AGREE * merged['listed_duration']
    return merged


def draw(workers, classes):
    fits = collect()
    ready = fits[fits['object_class'].isin([CLASSES[c] for c in classes]) & np.isfinite(fits['SNR']) & (fits['tau'] > 0)]
    todo = [(r.plot, int(r.tic), float(r.period), float(r.epoch), float(r.tau), float(r.SNR)) for r in ready.itertuples()
            if not os.path.exists(png_path(r.plot))]
    print(len(todo), 'figures to draw of', len(ready), 'fitted objects,', workers, 'workers', flush=True)
    t_start = time.time()
    for i, result in enumerate(pool_map(draw_one, todo, workers)):
        pd.DataFrame([result]).to_csv(place('log'), mode='a', header=not os.path.exists(place('log')), index=False)
        if result['error'] or (i + 1) % 200 == 0 or workers == 1:
            print('[%d/%d] %s %.0f s %s (%.1f min)' % (i + 1, len(todo), result['plot'], result['seconds'], result['error'], (time.time() - t_start) / 60), flush=True)
    print('done in %.1f min' % ((time.time() - t_start) / 60), flush=True)


def summary():
    fits = collect()
    fits['fitted'] = fits['SNR'].notna() if 'SNR' in fits else False
    fits['figure'] = [os.path.exists(png_path(p)) for p in fits['plot']]
    print(fits.groupby('object_class').agg(objects=('name', 'size'), fitted=('fitted', 'sum'), figure=('figure', 'sum')).to_string())
    if fits['fitted'].any():
        done = fits[fits['fitted']]
        print('fit stays at the listed ephemeris: %d of %d' % ((done['agrees'] == True).sum(), len(done)))
        print('SNR quantiles by review signal:')
        print(done.groupby(done['review_signal'].fillna('not reviewed'))['SNR'].quantile([.1, .5, .9]).unstack().round(1).to_string())
        print('errors:', fits['error'].fillna('').replace('', np.nan).dropna().value_counts().head(5).to_dict())


def install(apply):
    """Every planet and candidate of the targets, fitted or not, into the website's table of archive fits."""
    import site_figures
    fits = collect()
    planets = fits[fits['object_class'].isin([CLASSES['planet'], CLASSES['candidate']])].copy()
    planets['figure'] = [png_path(p) if os.path.exists(png_path(p)) else '' for p in planets['plot']]
    fitted = planets['SNR'].notna()
    print('planets and candidates %d: fitted %d (fit stays at the listed ephemeris %d), with a figure %d, errors %d'
          % (len(planets), fitted.sum(), (planets['agrees'] == True).sum(), (planets['figure'] != '').sum(),
             (planets['error'].fillna('') != '').sum()))
    if not apply:
        print('dry run; pass --apply to write known_fits.tsv and the figures')
        return
    planets.to_csv(place('site'), sep='\t', index=False)
    shown = planets['figure'] != ''
    if nea():
        shown &= planets['n_transits'] > 0                 # the catalog drops a planet without a transit in the data
    for r in planets[shown].itertuples():
        site_figures.install(r.figure, HERE, r.plot)
    print('installed', len(planets), 'rows; now run build_verified_catalogs.py')


if __name__ == '__main__':
    if '--nea' in sys.argv:
        sys.argv.remove('--nea')
        os.environ['KNOWN_FIGURES_GROUP'] = 'nea'
    stage = sys.argv[1] if len(sys.argv) > 1 else ''
    if stage == 'list':
        wanted = targets()
        print(wanted.groupby(['object_class', wanted['review_signal'].replace('', 'not reviewed')]).size().unstack(fill_value=0).to_string())
        print(len(wanted), 'objects on', wanted['tic'].nunique(), 'stars')
    elif stage == 'fit':
        fit(int(sys.argv[2]) if len(sys.argv) > 2 else 4, [int(a) for a in sys.argv[3:]] or None)
    elif stage == 'draw':
        draw(int(sys.argv[2]) if len(sys.argv) > 2 else 4, sys.argv[3:] or ['planet', 'candidate'])
    elif stage == 'summary':
        summary()
    elif stage == 'install':
        install('--apply' in sys.argv)
    else:
        sys.exit(__doc__)
