# Overview

We assign each candidate a significance by comparing its matched-filter SNR to a
per-star **null distribution**, built from null-simulation tests (NSTs). The
Gaussian p-value is

$$ p\text{-value} = 1 - \Phi\left(\frac{\mathrm{SNR} - \mu_\mathrm{null}}{\sigma_\mathrm{null}}\right), $$

where $\mu_\mathrm{null}$ and $\sigma_\mathrm{null}$ are the mean and standard
deviation of the maximum SNR recovered in each NST run on that star.

## How many NSTs are needed to estimate $\mu$ and $\sigma$?

Using ~100 NST runs for each of 100 randomly selected stars, we track the running
estimates of $\mu$ and $\sigma$ as a function of the number of NSTs used, $N$.

![Convergence of μ and σ with the number of NSTs](overview_nst_musigma.png)

## Is the Gaussian p-value reliable?

We pick, for each star, the SNR whose **true p-value is $0.01$** (the 99th
percentile of its null max-SNR samples) and estimate that p-value as a function
of $N$ from the Gaussian approximation.

![Gaussian p-value at the p ≈ 0.01 threshold](overview_nst_pvalue.png)

## How significant are the candidates?

Comparing the p-value distributions of the established (known) TOIs and our new
candidates — plotted as $-\log_{10}(p)$ on a log axis, so larger means more
significant:

![p-value distribution of known vs. new TOIs](overview_pvalue_dist.png)

## Habitable Zone

Each candidate's instellation flux is

$$ \frac{S}{S_\oplus} = \left(\frac{R_\star}{R_\odot}\right)^{2}\left(\frac{T_\mathrm{eff}}{5778~\mathrm{K}}\right)^{4} a^{-2}, $$

with $a$ [AU] $= M_\star^{1/3}(P/365.25)^{2/3}$ and $M_\star$ from $R_\star$ and
$\log g$. The conservative and optimistic habitable-zone boundaries are the
Kopparapu (2014) polynomial fits — runaway-greenhouse to maximum-greenhouse
(dark green) and recent-Venus to early-Mars (light green).

Markers are colored by $1 - p$-value (NST). Candidates with a dashed red ring
fail at least one systematics test (`spurious_transits`, `snrd`,
`num_transits`, or `sharp_freq`). The gold stars mark the three new candidates
sitting inside the optimistic HZ.

![Insolation–Teff diagram with HZ overlays](habitable_zone.png)

Of the **195** new candidates, **3** land inside the optimistic HZ
(TIC 30853470, 52307802, 178709444) and all 3 are also inside the conservative
HZ — but each fails at least one systematics test, so none currently pass.
For comparison, **23 of 5245** known TOIs are in the optimistic HZ (17 pass)
and **13** in the conservative HZ (11 pass).
