Known TESS planets and planet candidates that our search recovered on its own. Every row is a
signal from the TESS Objects of Interest (TOI) list which the pipeline found independently in
the `SPOC` light curves (it recovers $\gtrsim 99$ % of the TOIs at $\mathrm{SNR} > 20$ and
$95$ % of all of them). We keep signals with at least three transits, which the statistical
validation of [Robnik et al. (2024)](https://arxiv.org/abs/2407.17565) needs; that is about
80 % of the data, or about 5,000 stars.

**What changed in this version.** The list is now a verified one (see *Known planets and
binaries* in the Overview):

- **Class.** Each row is a *confirmed planet* or a *planet candidate*. The 724 signals that are
  eclipsing binaries, by the TOI disposition, the TESS-EB catalog or our own false-positive
  analysis, have moved to the Eclipsing Binaries page; 14 dispositioned as other false positives
  are no longer listed.
- **Signals fitted on the all-sky light curves.** 690 TOIs that the SPOC-based catalog did not
  hold were fitted at their verified ephemeris on the all-sky light curves (every TESS product
  of the star). They carry an SNR, a radius and a figure, but no p-value and no test results
  yet: the null searches of these stars are part of the all-sky run. They appear under
  *Show all catalog entries*.
- **Planets outside the TOI list.** The table opens on the TESS discoveries. Releasing the
  *TESS discoveries* button adds 2,991 confirmed planets that are in the NASA Exoplanet Archive
  but have no TOI number, mostly Kepler and K2 discoveries on stars that TESS also observed. They
  are listed at their archive ephemeris, which our search uses to mask them; we did not fit
  them, so they have no SNR, radius or figure.
- **Ephemeris.** Period, epoch and duration were checked against every TESS sector of the star.
  **Ephemeris from** says where the listed values come from: *fit and false-positive analysis*
  (our fitted period, with the epoch and duration refined by the posterior of the false-positive
  analysis), *review* or *single review* (read off the light curves by a reviewer, where our fit
  or another catalog had drifted), *catalog* or *fit* (no change). **Checked** says what the
  check found.
- **Duration** is the full width of the transit, from first to last contact. It is about 11 %
  longer than the centre-on-disc width that earlier versions of this table listed.
- **FPP**, **P(binary)** and **P(planet on another star)** are the probabilities of the
  false-positive analysis: FPP is everything except a planet on the target star. A candidate is
  listed here as long as its binary probability is below 50 %.

The SNR, the planet radius and the p-values come from our fit to the SPOC light curves. Where the
verified ephemeris had left an earlier fit by more than half a transit duration (416 rows), the
signal was fitted again at the verified ephemeris and its figure redrawn; the one row whose fit still
disagrees says so when clicked. Where the period moved to another period bin, the p-value uses the
null searches of the star's closest bin (8 rows, within a factor 2.5 in period); 3 rows whose period
changed by more have no p-value yet.

**Significance.** For each signal we report the $p$-value of its SNR against the star's own null
distribution in the signal's period bin (the null-search test of Robnik et al. 2024, with a
hierarchical prior pooled over all stars), multiplied by the star's number of period bins to
account for the search over all of them. This is the global $p$-value; the local one is listed
too (see the Overview).

**Tests.** By default the table shows only the signals that pass every automatic test:

- no transit better explained by a false-alarm scenario or by an overlapping deeper planet of the
  same star;
- a signal spread across its transits rather than carried by one;
- at least three transits;
- a folded dip shaped like a transit rather than like a sinusoid at the same period;
- an SNR that is significant against the star's own null distribution;
- no known eclipsing binary on the star.

**Show all catalog entries** reveals the rest, with the failing tests named. The tests are
deterministic and recomputed from the current fits; they are not a human judgement. Their
thresholds are calibrated on an injection and recovery run of the same pipeline: synthetic
planets of 2, 4 and 8 Earth radii are added to real light curves, and each threshold is the one
that recovers the most injected planets at a fixed false-alarm rate. Rows whose transit windows
contain no data, or whose transits are covered and show nothing, are not listed at all.
