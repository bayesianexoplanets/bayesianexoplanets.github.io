Transit timing variations (TTVs) and transit duration variations (TDVs) for Kepler candidates. All candidates from the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=cumulative), whose individual transits are significant enough (signal-to-noise ratio (SNR) > 3) and with period > 5 days are included. TTVs and TDVs were extracted as described in [Robnik et. al. 2026](https://academic.oup.com/mnras/article/547/3/stag419/8527731).

To identify significant TTVs we perform a $\chi^2$ test: we test the hypothesis that the signal-to-noise-ratio (SNR) of the planet's individual transits at their expected (periodic) emphemeris is consistent with no TTVs. Low p-value indicates that this is not the case and can be used to identify significant TTVs, we use p-value < 0.01 as a threshold.

Candidates that were identified by [Holczer et al. (2016)](https://iopscience.iop.org/article/10.3847/0067-0049/225/1/9) as having significant TTVs are also flagged.

| Column | Description |
|---|---|
| Kepler ID | Kepler Input Catalog identifier |
| KOI | Kepler Object of Interest number |
| Period [d] | Orbital period in days |
| τ [d] | Transit duration in days |
| SNR | Signal-to-noise ratio of the transit |
| p-value | SNRD test p-value for TTV significance |
| Holczer+2016 | Detected as TTV signal in Holczer et al. (2016) |
