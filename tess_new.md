New planet candidates found by our pipeline that are **not** in the TESS Objects of Interest
(TOI) catalog.

## Outlier detection: the most interesting candidates

We summarise each candidate by four features — planet radius, orbital period,
host $T_\mathrm{eff}$, and system multiplicity — and estimate the density of the
**known** TESS planet population in that space with a Gaussian KDE. Candidates
that land in low-density regions are flagged as out-of-distribution (OoD)
outliers (★) — the most interesting ones. The **Outlier score** column in the
table below ranks every candidate (higher = more unusual; flagged at an OoD
p-value $< 0.001$, $\approx 3\sigma$ — the most extreme 0.1% of the
known-population density). 19 candidates are flagged interesting.

![New candidates vs. known TESS planets in feature space](new_outliers.png)

Most of the outliers are ultra-short-period ($P \lesssim 0.3$ d) or long-period
($P \gtrsim 180$ d) planets, unusually large/small radii, or members of
multi-planet systems (e.g. TIC 178155732 = HR 858).

| TIC | Period [d] | Radius [R⊕] | Teff [K] | Mult. | Outlier score |
|---|---|---|---|---|---|
| 420112589 | 282.83 | 1.57 | 3212 | 1 | 7.35 |
| 231077395 | 0.234 | 0.67 | 6339 | 1 | 7.10 |
| 232635922 | 183.21 | 2.48 | 3601 | 1 | 4.89 |
| 90127880 | 0.247 | 5.32 | 8196 | 1 | 4.72 |
| 278862750 | 0.260 | 1.06 | 6775 | 1 | 4.53 |
| 280206394 | 0.296 | 1.15 | 6236 | 1 | 3.78 |
| 395171208 | 0.231 | 1.47 | 6839 | 1 | 3.77 |
| 326386495 | 0.227 | 0.73 | 4781 | 1 | 3.59 |
| 30853470 | 233.84 | 2.22 | 4613 | 1 | 3.44 |
| 161477033 | 0.813 | 1.91 | 5245 | 2 | 3.26 |
| 178155732 | 20.49 | 1.64 | 6309 | 3 | 3.12 |
| 169191182 | 0.311 | 1.43 | 6180 | 1 | 2.96 |
| 52928939 | 0.217 | 1.38 | 5417 | 1 | 2.87 |
| 112115898 | 1.310 | 0.54 | 3190 | 1 | 2.86 |
| 287328202 | 0.295 | 1.50 | 6088 | 1 | 2.82 | 
