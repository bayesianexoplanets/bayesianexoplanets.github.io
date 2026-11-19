"""The two verified TESS catalogs of the website, from the pipeline's own catalog (`tois.csv`) and the verified
known-object table of the all-sky search (`analysis/known_objects.py verify`; logs/DESIGN.md 2026-11-16 to 2026-11-19).

    python build_verified_catalogs.py

- `tess_planets.csv`: the rows of `tois.csv` whose object is a planet or planet candidate, with the verified period,
  epoch and full duration; the fit's own values stay in `Period_fit`, `Epoch_fit`, `Duration_fit`, `Phase_fit`,
  `Tau_fit`. `fit_stale` marks rows whose SNR, radius, p-value and figure were computed at an ephemeris that the
  verified one leaves by more than half a fitted duration inside the star's observed sectors.
- `tess_ebs.csv`: every verified eclipsing binary (TESS-EB catalog, TOIs dispositioned as EBs, candidates reclassified
  by the false-positive analysis or a review), with its secondary eclipse where one is known.
Rows refitted at their verified ephemeris (`refit_verified.py`, `refits.tsv`) take the refit's SNR, radius, transit count
and fit columns, and the p-value of the new SNR in the stored null distribution of their period bin; a row whose period
left its bin borrows the star's stored null closest in period (`null_bin_period`), up to a factor CLOSEST_BIN.
A second single-agent look at ten refitted rows (`known_review_refit/findings.tsv`) sets their ephemeris the same way,
except where its verdict is `unclear`; `planet_notes.tsv` holds the one-line remark shown with a row (`Note`).
`tess_planets.csv` also holds two kinds of rows that the website's own pipeline never fitted:
- the verified planets and candidates fitted on the all-sky archive (`known_figures.py install` -> `known_fits.tsv`;
  `Data` = `archive`): SNR, radius and figure from that fit, no null searches and therefore no p-value and no tests;
- the confirmed planets without a TOI number (`NEA_only`): table objects that hold a NASA Exoplanet Archive planet and neither a
  TOI nor a TESS-EB entry, on a star of the archive. The search masks them at their catalog values without fitting
  them, so they carry the catalog ephemeris only; the website hides them until "TESS discoveries" is switched off.
`plot` is the figure of a row (the per-star index of `tois.csv` is fixed here, so removing rows cannot shift it).
Text cells never contain a comma: the website splits lines on commas. `site_stats.json` holds the counts the pages quote.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from analysis import known_objects as ko
import sm_pvalue
import site_figures

TABLE = ko.OUT_ROOT + 'known_objects_20261114_verified_v2.tsv'
FINDINGS = ko.scratch + 'allsky_search/known_review_opus/findings.tsv'   # the single-agent check of 80 candidates
REFIT_FINDINGS = ko.scratch + 'allsky_search/known_review_refit/findings.tsv'   # the second look at ten refitted rows
NOTES = os.path.join(HERE, 'planet_notes.tsv')
KNOWN_FITS = os.path.join(HERE, 'known_fits.tsv')                       # known_figures.py install
EB_FIGURES = os.path.join(HERE, 'eb_figures.tsv')                       # eb_figure.py install
EB_REVIEW = os.path.join(HERE, 'eb_review2.tsv')                        # allsky_search/known_review_eb/apply_review.py --write
SECOND_REVIEW = {'clear': 'second review: clear signal', 'weak': 'second review: weak signal', 'none': 'second review: no signal'}
PLANET_CLASSES = {'planet': 'confirmed planet', 'planet candidate': 'planet candidate'}
STALE = 0.5                                                             # fitted durations (ko.DISPUTE)
REVIEW = {'clear': 'review: clear signal', 'weak': 'review: weak signal', 'none': 'review: no signal'}
REFITS = os.path.join(HERE, 'refits.tsv')
SAME_BIN = 0.01                                                         # relative period change that keeps the null's period bin
CLOSEST_BIN = 2.5                                                       # largest period factor to a borrowed null (user, 2026-11-19)
NULL_COLUMNS = ['sm_sf_grid', 'nst_samples', 'μ(SNR | null)', 'σ(SNR | null)', 'n_period_bins']


class VerifiedCatalogs:
    """The website rows matched to the verified table, and the two catalogs built from them."""

    def __init__(self):
        from TESS.archive_census import ArchiveCensus
        self.site = pd.read_csv(os.path.join(HERE, 'tois.csv'), low_memory=False)
        self.table = pd.read_csv(TABLE, sep='\t', low_memory=False)
        for column in ('review', 'review_signal', 'posterior'):
            self.table[column] = self.table[column].fillna('')
        second = pd.read_csv(REFIT_FINDINGS, sep='\t')
        second = second[second['verdict'] != 'unclear'].assign(signal=lambda f: f['signal_archive'])   # every sector, as the other checks
        self.findings = pd.concat([pd.read_csv(FINDINGS, sep='\t'), second], sort=False).set_index('candidate')
        assert self.findings.index.is_unique, 'a candidate has two single reviews'
        self.notes = pd.read_csv(NOTES, sep='\t').set_index('candidate')['note'].to_dict() if os.path.exists(NOTES) else {}
        self.fpp = pd.read_csv(ko.POSTERIORS + 'fpp_summary.csv').set_index('candidate')
        self.families = ko.scenario_families()
        self.windows = ko.ObservedWindows(ArchiveCensus.load(ko.CENSUS_PATH, ko.ARCHIVE_ROOT))
        self.refits = (pd.read_csv(REFITS, sep='\t').set_index('plot') if os.path.exists(REFITS) else pd.DataFrame())
        self.owner = {}
        for index, tic, members in zip(self.table.index, self.table['tic'], self.table['members'].astype(str)):
            for name in members.split('|'):
                self.owner.setdefault((int(tic), name), index)

        site = self.site
        site['plot'] = 'plots/' + site['TIC'].astype(str) + '/' + site.groupby('TIC').cumcount().astype(str) + '.jpg'
        name = np.where(site['TOI'].notna(), 'TOI ' + site['TOI'].astype(str), 'TIC ' + site['TIC'].astype(str) + ' website')
        site['candidate'] = np.where(site['TOI'].notna(), 'TOI-' + site['TOI'].astype(str), '')
        site['row'] = [self.owner.get((int(tic), n), -1) for tic, n in zip(site['TIC'], name)]


    def checked(self, row):
        """How the table object was checked, in words."""
        if row['review_signal'] in REVIEW:
            return REVIEW[row['review_signal']]
        return 'false-positive analysis' if row['posterior'] != '' else 'not checked'


    def source(self, row):
        """Where the table object's ephemeris comes from, in words."""
        return ('fit and false-positive analysis' if row['posterior'] != '' else
                'review' if row['review_signal'] in ko.REVIEW_SIGNAL else 'catalog')


    def planet_rows(self):
        """`tess_planets.csv`: one row per website row of planet class, at its verified ephemeris."""
        rows = []
        for _, s in self.site.iterrows():
            if s['row'] < 0:
                continue
            t = self.table.loc[s['row']]
            finding = self.findings.loc[s['candidate']] if s['candidate'] in self.findings.index else None
            own = abs(s['Period'] / t['period'] - 1.) < ko.SAME_PERIOD         # the table object is this row's signal
            object_class = t['object_class']
            if finding is not None and finding['verdict'] == 'both_real_two_signals':
                object_class = finding['class']                                # the object holds the other signal
            if object_class not in PLANET_CLASSES:
                continue
            if finding is not None:
                period, epoch, duration = float(finding['period']), float(finding['epoch']), float(finding['duration_h']) / 24.
                checked, source = 'single review: %s signal' % finding['signal'], 'single review'
            elif own:
                period, epoch, duration = float(t['period']), float(t['epoch']), float(t['duration'])
                checked, source = self.checked(t), self.source(t)
            else:
                period, epoch, duration = float(s['Period']), float(s['Epoch']), float(s['Duration'])
                checked, source = 'not checked', 'fit'
            gap = (ko.transit_gap(dict(period=s['Period'], epoch=s['Epoch']), dict(period=period, epoch=epoch), self.windows.of(int(s['TIC'])))
                   if abs(period / s['Period'] - 1.) < ko.SAME_PERIOD else np.inf)
            row = s.drop(labels=['row', 'candidate']).to_dict()
            fitted = source in ('fit', 'fit and false-positive analysis')
            row.update(Period=period, Epoch=epoch, Duration=duration, Class=PLANET_CLASSES[object_class], Checked=checked,
                       Ephemeris=source, fit_stale=bool(gap > STALE * s['Duration']), null_bin_period=np.nan,
                       Period_fit=s['Period'], Epoch_fit=s['Epoch'], Duration_fit=s['Duration'], Phase_fit=s['Phase'], Tau_fit=s['Tau'],
                       err_Period=s['err_Period'] if fitted else np.nan, err_Epoch=s['err_Epoch'] if fitted else np.nan,
                       err_Duration=np.nan)
            in_analysis = s['candidate'] in self.fpp.index and (own or finding is not None)
            row.update(FPP=self.fpp.at[s['candidate'], 'FPP'] if in_analysis else np.nan,
                       P_binary=self.families.at[s['candidate'], 'binary'] if in_analysis else np.nan,
                       P_elsewhere=self.families.at[s['candidate'], 'elsewhere'] if in_analysis else np.nan,
                       Note=self.notes.get(s['candidate'], ''), Data='SPOC', NEA_only=False)
            del row['Phase'], row['Tau']
            if row['fit_stale'] and s['plot'] in self.refits.index and self.refits.at[s['plot'], 'agrees'] == True:
                row.update(self.refit(row, self.refits.loc[s['plot']]))
            rows.append(row)
        return pd.DataFrame(rows)


    def archive_rows(self, listed):
        """The verified planets and candidates without a website fit, as fitted on the all-sky archive.

        `listed` holds the (TIC, TOI) pairs that already have a row: a TOI is never listed twice, and an object without
        a TOI number (the website's own signal at another period) is left to its website row.
        """
        if not os.path.exists(KNOWN_FITS):
            return pd.DataFrame()
        fits = pd.read_csv(KNOWN_FITS, sep='\t')
        table = self.table.set_index(['tic', 'name'])
        rows = []
        for f in fits.itertuples():
            t = table.loc[(f.tic, f.name)]
            toi = f.toi if np.isfinite(f.toi) else member_toi(t)
            if not np.isfinite(toi) or (f.tic, toi) in listed:
                continue
            fitted = np.isfinite(f.SNR)
            rows.append({
                'TIC': f.tic, 'TOI': toi, 'known_as': '', 'Class': PLANET_CLASSES[t['object_class']],
                'Period': t['period'], 'Epoch': t['epoch'], 'Duration': t['duration'], 'Checked': self.checked(t), 'Ephemeris': self.source(t),
                'SNR': f.SNR, 'Radius_planet': f.radius, 'Radius_planet_errp': f.radius_errp, 'Radius_planet_errm': f.radius_errm,
                'Radius': f.star_radius, 'Number of Valid Transits': f.n_transits if fitted and f.n_transits >= 0 else np.nan,
                'Period_fit': f.fit_period, 'Epoch_fit': f.fit_epoch, 'Duration_fit': 2. * f.tau, 'Tau_fit': f.tau, 'Phase_fit': f.phase,
                'fit_stale': bool(fitted and f.agrees != True), 'plot': f.plot if isinstance(f.figure, str) and f.figure != '' else '',
                'Note': '', 'Data': 'archive', 'NEA_only': False})
        return pd.DataFrame(rows)


    def nea_rows(self):
        """The confirmed planets without a TOI number on archive stars: catalog values, never fitted."""
        t = self.table
        shown = set(self.site.loc[self.site['row'] >= 0, 'row'])
        nea = t[t['frozen'] & t['masked'] & t['in_archive'] & t['object_class'].isin(PLANET_CLASSES) & ~t.index.isin(shown)]
        return pd.DataFrame({
            'TIC': nea['tic'].to_numpy(), 'known_as': [planet_name(r) for _, r in nea.iterrows()],
            'Class': nea['object_class'].map(PLANET_CLASSES).to_numpy(), 'Period': nea['period'].to_numpy(), 'Epoch': nea['epoch'].to_numpy(),
            'Duration': nea['duration'].to_numpy(), 'Checked': [self.checked(r) for _, r in nea.iterrows()], 'Ephemeris': 'catalog',
            'plot': '', 'Note': '', 'Data': '', 'NEA_only': True})


    def refit(self, row, fit):
        """The columns a refit at the listed ephemeris replaces (as `merge_known_run.py` does for a known-planet run)."""
        new = dict(SNR=fit['SNR'], Period_fit=fit['period'], Epoch_fit=fit['epoch'], Duration_fit=2. * fit['tau'], Tau_fit=fit['tau'],
                   Phase_fit=fit['phase'], fit_stale=False)
        if np.isfinite(fit['radius']):
            new.update(Radius_planet=fit['radius'], Radius_planet_errp=fit['radius_errp'], Radius_planet_errm=fit['radius_errm'])
        if fit['n_transits'] >= 0:
            new['Number of Valid Transits'] = int(fit['n_transits'])
        if abs(row['Period'] / row['Period_fit'] - 1.) >= SAME_BIN:           # the period left its bin: the star's closest stored null
            star = self.site[(self.site['TIC'] == row['TIC']) & self.site['sm_sf_grid'].notna()]
            distance = np.abs(np.log(row['Period'] / star['Period']))
            if len(star) and np.exp(distance.min()) <= CLOSEST_BIN:
                source = star.loc[distance.idxmin()]
                new.update({column: source[column] for column in NULL_COLUMNS})
                new['null_bin_period'] = source['Period']
            else:
                new.update({column: np.nan for column in NULL_COLUMNS[:4]})      # no null search near this period yet
        grid = new.get('sm_sf_grid', row['sm_sf_grid'])
        bins = new.get('n_period_bins', row['n_period_bins'])
        local = sm_pvalue.log10p_from_grid(grid, fit['SNR']) if isinstance(grid, str) else np.nan
        new['log10(p local)'] = local
        new['log10(p value)'] = min(0., local + np.log10(bins)) if np.isfinite(local) and bins > 0 else np.nan
        return new


    def eb_rows(self):
        """`tess_ebs.csv`: one row per verified eclipsing binary of the table."""
        t = self.table
        ebs = t[(t['object_class'] == 'EB') & t['masked']].copy()
        site = self.site[self.site['row'] >= 0].drop_duplicates('row').set_index('row')
        on_site = site.reindex(ebs.index)
        own = (on_site['Period'] / ebs['period'] - 1.).abs() < ko.SAME_PERIOD
        members = ebs['members'].astype(str).str.split('|')
        toi = members.apply(lambda m: next((x[4:] for x in m if x.startswith('TOI ') and x[4:].replace('.', '').isdigit()), ''))
        catalog = members.apply(lambda m: next((x[8:] for x in m if x.startswith('TESS-EB ')), ''))
        origin = np.where(ebs['object_class_catalog'] == 'EB',
                          np.where(catalog != '', 'TESS-EB catalog', 'TOI disposition'),
                          np.where(ebs['fpp_binary'] > ko.FPP_MAX, 'false-positive analysis', 'review'))
        name = np.where(toi != '', 'TOI-' + toi, np.where(catalog != '', 'TESS-EB ' + catalog, ebs['name'].astype(str).str.replace(',', ' ')))
        return pd.DataFrame({
            'TIC': ebs['tic'].to_numpy(), 'Name': name, 'TOI': toi.to_numpy(), 'Catalog': catalog.to_numpy(), 'Origin': origin,
            'Period': ebs['period'].to_numpy(), 'Epoch': ebs['epoch'].to_numpy(), 'Duration': ebs['duration'].to_numpy(),
            'Secondary_phase': ebs['secondary_phase'].to_numpy(), 'Secondary_duration': ebs['secondary_duration'].to_numpy(),
            'Checked': [self.checked(r) for _, r in ebs.iterrows()], 'P_binary': ebs['fpp_binary'].to_numpy(),
            'SNR': np.where(own, on_site['SNR'], np.nan), 'in_archive': ebs['in_archive'].to_numpy(),
            'plot': np.where(own, on_site['plot'], '')})


    def second_review(self, ebs):
        """The second review of the binaries (on the EB figures' detrending): its period, epoch, width, secondary eclipse
        and verdict replace the row's. Without a secondary of its own, the listed one is kept only where the review
        left the period and epoch as listed."""
        if not os.path.exists(EB_REVIEW):
            return ebs
        review = pd.read_csv(EB_REVIEW, sep='\t').drop_duplicates(['Name', 'TIC']).set_index(['Name', 'TIC'])
        key = pd.MultiIndex.from_frame(ebs[['Name', 'TIC']])
        hit = key.isin(review.index)
        found = review.loc[key[hit]]
        ebs = ebs.copy()
        own_secondary = found['secondary_phase'].notna().to_numpy()
        unmoved = found['unmoved'].to_numpy(dtype=bool)
        for column, source in (('Period', 'period'), ('Epoch', 'epoch'), ('Duration', 'duration')):
            ebs.loc[hit, column] = found[source].to_numpy()
        for column, source in (('Secondary_phase', 'secondary_phase'), ('Secondary_duration', 'secondary_duration')):
            ebs.loc[hit, column] = np.where(own_secondary, found[source].to_numpy(), np.where(unmoved, ebs.loc[hit, column].to_numpy(), np.nan))
        ebs.loc[hit, 'Checked'] = np.where(found['verdict'] == 'not_eclipse', 'second review: not an eclipse', found['signal'].map(SECOND_REVIEW))
        return ebs


    def eb_figures(self, ebs):
        """Each binary's EB figure (`eb_figure.py`) where it is installed, with the kind of curve it shows."""
        figure = ('plots_eb/' + ebs['TIC'].astype(str) + '/' + ebs.groupby('TIC').cumcount().astype(str) + '.webp').to_numpy()
        there = np.array([os.path.exists(os.path.join(HERE, f)) for f in figure])
        drawn = (pd.read_csv(EB_FIGURES, sep='\t').set_index('figure') if os.path.exists(EB_FIGURES)
                 else pd.DataFrame(columns=['route', 'period_orb'], dtype=float))
        curve = pd.Series(figure).map(drawn['route']).fillna('').to_numpy()
        curve = np.where((curve == 'none') & ~ebs['Checked'].str.endswith('no signal').to_numpy(), 'no curve', curve)   # a fit that missed
        return ebs.assign(plot=np.where(there, figure, ebs['plot']),
                          Curve=np.where(there, curve, 'transit fit'),
                          Fold_period=np.where(there, pd.Series(figure).map(drawn['period_orb']).to_numpy(dtype=float), np.nan))


def member_toi(row):
    """The TOI number among the object's members, NaN without one."""
    numbers = [m[4:] for m in str(row['members']).split('|') if m.startswith('TOI ') and m[4:].replace('.', '').isdigit()]
    return float(numbers[0]) if numbers else np.nan


def planet_name(row):
    """The planet's catalog name: the first member that carries a reference tag, without the tag."""
    tagged = [m for m in str(row['members']).split('|') if '[' in m]
    return (tagged[0].split('[')[0] if tagged else str(row.name[1] if isinstance(row.name, tuple) else row['name'])).strip().replace(',', ' ')


def write(frame, name):
    frame = frame.assign(plot=[site_figures.published(HERE, p) or (p if isinstance(p, str) else '') for p in frame['plot']])
    text = frame.select_dtypes(include='object')
    assert not text.apply(lambda c: c.astype(str).str.contains(',').any()).any(), 'a text cell contains a comma'
    frame.to_csv(os.path.join(HERE, name + '.tmp'), index=False)
    os.replace(os.path.join(HERE, name + '.tmp'), os.path.join(HERE, name))
    print('%s: %d rows' % (name, len(frame)))


if __name__ == '__main__':
    catalogs = VerifiedCatalogs()
    planets, ebs = catalogs.planet_rows(), catalogs.eb_rows()
    ebs = catalogs.eb_figures(catalogs.second_review(ebs))
    ebs.loc[ebs['plot'] == '', 'Curve'] = ''
    archive, nea = catalogs.archive_rows(set(zip(planets['TIC'], planets['TOI']))), catalogs.nea_rows()
    write(pd.concat([planets, archive, nea], ignore_index=True, sort=False)[list(planets.columns)], 'tess_planets.csv')
    print('rows fitted on the SPOC files %d, on the archive %d (with SNR %d, with a figure %d), confirmed planets without a TOI number %d on %d stars'
          % (len(planets), len(archive), archive['SNR'].notna().sum() if len(archive) else 0, (archive['plot'] != '').sum() if len(archive) else 0,
             len(nea), nea['TIC'].nunique()))
    write(ebs, 'tess_ebs.csv')
    print('planets: class %s; passing %d; stale fit %d (passing %d); ephemeris %s; checked %s'
          % (planets['Class'].value_counts().to_dict(), planets['passed_all_tests'].sum(), planets['fit_stale'].sum(),
             (planets['fit_stale'] & planets['passed_all_tests']).sum(), planets['Ephemeris'].value_counts().to_dict(),
             planets['Checked'].value_counts().to_dict()))
    print('EBs: origin %s; checked %s; with a secondary %d; with a figure %d; in the archive %d'
          % (ebs['Origin'].value_counts().to_dict(), ebs['Checked'].value_counts().to_dict(), ebs['Secondary_phase'].notna().sum(),
             (ebs['plot'] != '').sum(), ebs['in_archive'].sum()))
    new = pd.read_csv(os.path.join(HERE, 'tois_new.csv'), low_memory=False)
    signal = lambda frame: int(frame['Checked'].str.contains('clear|weak').sum())
    listed = pd.concat([planets, archive], ignore_index=True, sort=False)      # the TESS discoveries: SPOC and archive fits
    stats = dict(updated=time.strftime('%Y-%m-%d'),
                 planets=len(listed), confirmed=int((listed['Class'] == 'confirmed planet').sum()),
                 candidates=int((listed['Class'] == 'planet candidate').sum()), planets_passing=int(planets['passed_all_tests'].sum()),
                 planets_spoc=len(planets),
                 planets_reviewed_signal=signal(planets), planets_analysis_only=int((planets['Checked'] == 'false-positive analysis').sum()),
                 planets_stale=int(planets['fit_stale'].sum()), planet_hosts=int(listed['TIC'].nunique()),
                 planets_refitted=int(planets['plot'].isin(catalogs.refits.index[catalogs.refits['agrees'] == True]).sum()) if len(catalogs.refits) else 0,
                 ebs=len(ebs), ebs_catalog=int((ebs['Origin'] == 'TESS-EB catalog').sum()), ebs_toi=int((ebs['Origin'] == 'TOI disposition').sum()),
                 ebs_reclassified=int(ebs['Origin'].isin(['false-positive analysis', 'review']).sum()), ebs_reviewed_signal=signal(ebs),
                 ebs_no_signal=int(ebs['Checked'].str.endswith('no signal').sum()),
                 ebs_second_review=int(ebs['Checked'].str.startswith('second review').sum()), ebs_secondary=int(ebs['Secondary_phase'].notna().sum()),
                 ebs_figures=int((ebs['plot'] != '').sum()), eb_hosts=int(ebs['TIC'].nunique()),
                 new_candidates=len(new), new_passing=int(new['passed_all_tests'].sum()),
                 planets_archive=len(archive), planets_other_surveys=len(nea), other_survey_hosts=int(nea['TIC'].nunique()),
                 ebs_curves=int(ebs['Curve'].isin(['posterior', 'physical fit', 'two-Gaussian fit']).sum()),
                 new_figure_ext='.webp' if all(os.path.exists(os.path.join(HERE, 'plots_new/%d/%d.webp' % (t, i)))
                                               for t, i in zip(new['TIC'], new['_cand_idx'])) else '.jpg')
    with open(os.path.join(HERE, 'site_stats.json'), 'w') as f:
        json.dump(stats, f, indent=1)
    print(json.dumps(stats))
