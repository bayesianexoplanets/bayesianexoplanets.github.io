# Detection Pipeline Overview

This page summarises the automated pipeline that searches TESS (and Kepler) light
curves for transiting exoplanets. The full methodology is described in
[Robnik, Seljak, Jenkins & Bryson (2026)](https://arxiv.org/abs/2601.07465).

---

## Why is this hard?

Detecting Earth-sized planets in habitable-zone orbits requires finding
transit depths of order **100 ppm** (0.01 % of stellar flux) against a
background of stellar variability that is typically 2–5× larger. A single
transit event lasts only a few hours in an orbital period of ~200–400 days,
so a multi-year baseline yields only a handful of events per star. Any
instrumental artefact — a thermal transient, a pixel sensitivity drop, a
cosmic ray — can mimic a genuine transit. The pipeline must therefore be both
sensitive enough to detect real signals near the noise floor *and* robust
enough to reject the inevitable flood of false alarms.

---

## Data model

The fundamental assumption is that the observed flux decomposes as

$$
F(t) = F_\mathrm{GP}(t) + \epsilon(t) + \Delta(t),
$$

where $F_\mathrm{GP}$ is the **stellar variability** (a correlated Gaussian
process), $\epsilon$ is uncorrelated Gaussian noise whose variance changes
with the local flux, and $\Delta$ is the planet transit signal. Working in
Fourier space, the stellar variability is fully characterised by its **power
spectral density** $P_k$. Given $P_k$ the data can be **whitened**: after
whitening the planet signal is a template whose shape is determined by the
transit physics alone, and the noise is approximately white.

The power spectrum $P_k$ is estimated directly from the data by smoothing the
periodogram with an adaptive window whose width grows with frequency to
account for possible non-stationarity. Because stellar variability appears as
a broad, smooth feature in the power spectrum, this estimate is stable even
in the presence of a weak planet transit.

---

## Pipeline overview

The pipeline runs in three sequential stages:

![Pipeline stages: preprocessing → detection → vetting](pipeline_overview.png)

---

## Stage 1 — Preprocessing

Before searching for planets, several sources of non-Gaussianity must be
removed.

### Outlier Gaussianization

Isolated outliers — bright pixels, saturation spikes — would dominate the
matched filter if left in place. The **Gaussianization** step applies a
monotonic transformation that maps the empirical flux distribution to a
standard Gaussian while preserving extended correlations (such as a planet
transit dip). The transformation is iterated: estimate the power spectrum,
apply the transformation, repeat. This is done in 2–3 passes.

### Known-planet masking and TTV fitting

If a star has previously identified large planets, their transit signals
contaminate the power-spectrum estimate and suppress sensitivity to smaller
companions. For each known planet, individual transits are tested for
significance; those with SNR > 8 are refitted with a **transit timing
variation (TTV) model** to measure timing offsets. All known transits are
then masked — replaced by the best-fit Gaussian-process model — before the
subsequent search.

### Localized defect removal

Instrumental events such as sudden pixel sensitivity drops and thermal
transients produce localised discontinuities in the light curve. A
**template bank** of step-function and exponential-decay models is evaluated
at every cadence; any site with an SNR > 8 match is flagged and replaced by
the GP interpolant. Without this step, a single defect can generate many
harmonic peaks in the detection grid.

### Harmonic / notch filtering

Many stars show narrow periodic peaks in their power spectra — from
spacecraft electronics, reaction-wheel harmonics, or pulsation overtones.
These are identified via wavelet analysis and **notch-filtered**: the
matched-filter sensitivity at their frequencies is suppressed before the
period grid search begins.

---

## Stage 2 — Detection

### Matched-filter test statistic

For a transit with period $P$, epoch $\phi$, and duration $\tau$, the
**matched-filter SNR** is

$$
\mathrm{SNR}(P,\phi,\tau) = \frac{\sum_k w_k \, \delta_k}{\sqrt{\sum_k w_k}},
$$

where $\delta_k$ is the whitened flux in the $k$-th cadence and $w_k$ is the
matched-filter weight (proportional to the transit template). The optimal
Bayesian test statistic is the **Bayes factor** — the ratio of evidence with
a planet to evidence without one — which integrates over the transit
amplitude. For Gaussian noise its logarithm is approximately
$\frac{1}{2}\mathrm{SNR}^2$ minus a depth-prior penalty.

### Period–phase–duration grid

The Bayes factor is evaluated on a dense grid:
- **Period:** 0.2–40 days with ~150 000 trial periods spaced so that the
  phase error accumulated over the baseline never exceeds half a cadence.
- **Phase:** For each period, the phase that maximises the Bayes factor is
  found analytically in $O(N\log N)$ time using a Fourier-space inner product.
- **Duration:** The transit duration $\tau$ is constrained by Kepler's third
  law given the stellar density $\rho_\star$, with a ±15 % uncertainty folded
  into a prior. Only durations consistent with a physical orbit are evaluated.

The grid search identifies **Threshold Crossing Peaks (TCPs)** — local maxima
above a detection threshold. An iterative masking step converts TCPs to
independent **Threshold Crossing Events (TCEs)**: the strongest peak is
accepted, its signal is subtracted, and the grid is re-searched until no
further peaks exceed the threshold.

---

## Stage 3 — Vetting

Reaching the noise floor means many TCEs are instrumental false alarms. Three
complementary tests reject them.

### Individual-transit glitch test

Each transit event is tested against a **localized-defect hypothesis**: is
the flux dip better explained by an instrumental artefact than by a transit
template? If the false-alarm model fits significantly better (empirical
p-value < 0.1) the transit is flagged *spurious*. If more than 50 % of the
total $\mathrm{SNR}^2$ is carried by spurious transits, the entire TCE is
rejected.

### Gap-proximity test

Transits within ~0.3 days of a data gap are flagged as potentially spurious.
Data gaps in TESS often coincide with thermal settling events, so cadences
near a gap boundary are less reliable. The same 50 % threshold applies.

### Per-transit SNR consistency (χ² test)

For a genuine planetary transit, each event should contribute to the total
detection with an SNR proportional to its flux weight. The pipeline tests
whether the distribution of per-transit SNR values is consistent with this
expectation using a **χ² uniformity test**. A very non-uniform distribution —
one or two anomalously deep events driving the entire signal — flags an
unmodeled systematic. TCEs with p-value < 0.01 are rejected.

---

## Significance: null-simulation tests (NSTs)

The matched-filter SNR depends on each star's noise properties in a
complicated way, making it difficult to assign a universal detection threshold.
Instead, the pipeline runs **null-simulation tests (NSTs)**: ten independent
searches per star with a slightly perturbed search grid (different random seed),
so the pipeline finds noise peaks rather than real planets. The **empirical
null distribution** of max-SNR over these runs defines the per-star noise floor.

The significance of a candidate at $\mathrm{SNR}_\mathrm{cand}$ is

$$
p = 1 - \Phi\!\left(\frac{\mathrm{SNR}_\mathrm{cand} - \mu_\mathrm{null}}{\sigma_\mathrm{null}}\right),
$$

where $\mu_\mathrm{null}$ and $\sigma_\mathrm{null}$ are the mean and standard
deviation of the NST null distribution.

**How many NSTs are needed?** Using ~100 NST runs per star we track the
convergence of $\mu$ and $\sigma$ as a function of $N$:

![Convergence of μ and σ with the number of NSTs](overview_nst_musigma.png)

Both estimates stabilise within ~5 runs; the ten used in production give a
reliable Gaussian approximation. The NST-based p-value at the $p \approx 0.01$
significance level converges to the correct value within 5–8 runs:

![Gaussian p-value at the p ≈ 0.01 threshold vs number of NSTs](overview_nst_pvalue.png)

**How significant are the candidates?**

Comparing the $-\log_{10}(p)$ distributions of known TOIs and our new candidates
(larger = more significant):

![p-value distribution of known TOIs vs new candidates](overview_pvalue_dist.png)

The new candidates are systematically less significant than confirmed TOIs —
as expected near the detection limit — but their p-values are well below the
per-star noise floor.

---

## Habitable zones

Each candidate's **insolation flux** relative to Earth is

$$
\frac{S}{S_\oplus} = \left(\frac{R_\star}{R_\odot}\right)^{2}
\left(\frac{T_\mathrm{eff}}{5778\,\mathrm{K}}\right)^{4} a^{-2},
$$

with semi-major axis $a\,[\mathrm{AU}] = M_\star^{1/3}(P/365.25)^{2/3}$ and
$M_\star$ inferred from $R_\star$ and $\log g$. The habitable-zone boundaries
follow [Kopparapu et al. (2014)](https://doi.org/10.1088/2041-8205/787/2/L29):
*runaway greenhouse* to *maximum greenhouse* (dark green, conservative HZ) and
*recent Venus* to *early Mars* (light green, optimistic HZ).

### TESS new candidates

Markers are coloured by $\log_{10}(p\text{-value})$, clipped to $[-3, 0]$.
A dashed red ring marks candidates that fail at least one vetting test
(`spurious`, `snrd`, or `ntransits`).
Gold stars mark candidates inside the optimistic HZ.

![TESS candidates: insolation–Teff diagram with HZ overlays](habitable_zone.png)

Of the **195** new candidates, **3** land inside the optimistic HZ
(TIC 30853470, 52307802, 178709444) — all also inside the conservative HZ —
but each fails at least one vetting test, so none currently pass.
For comparison, 23 of 5 245 known TOIs are in the optimistic HZ (17 pass) and
13 in the conservative HZ (11 pass).

### Kepler candidates (Robnik et al. 2026)

The same pipeline applied to the Kepler dataset:

![Kepler candidates: insolation–Teff diagram with HZ overlays](kepler_hz.png)

Several long-period habitable-zone candidates found by the pipeline (red circles)
were not in previously published candidate lists (blue, eliminated by earlier
vetting). KOI 8063.01, 8107.01, and 8242.01 are identified as false alarms
by the refined vetting stage.
