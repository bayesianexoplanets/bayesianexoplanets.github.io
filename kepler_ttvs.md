# Kepler TTVs

Transit timing variations (TTVs) detected in Kepler light curves using a Bayesian SNRD test.
The catalog contains 1,016 planet candidates. The **p-value** column gives the probability of observing
the measured SNRD under the null hypothesis of no TTVs; small values indicate significant TTV signal.
**Holczer+2016** marks candidates also reported in the [Holczer et al. (2016)](https://doi.org/10.3847/0067-0049/225/1/9) TTV catalog.

| Column | Description |
|---|---|
| Kepler ID | Kepler Input Catalog identifier |
| KOI | Kepler Object of Interest number |
| Period [d] | Orbital period in days |
| Phase [d] | Reference transit time (days) |
| τ [d] | Transit duration in days |
| SNR | Signal-to-noise ratio of the transit |
| p-value | SNRD test p-value for TTV significance |
| Holczer+2016 | Detected as TTV signal in Holczer et al. (2016) |
