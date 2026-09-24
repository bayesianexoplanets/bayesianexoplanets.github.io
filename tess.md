Here we compiled a list of all exoplanet candidates, independently verified using our
pipeline (we recovered $\gtrsim 99$ % at $\mathrm{SNR} > 20$ and $95$ % of all candidates),
which were already known from the TESS Objects of Interest (TOI) catalog. For consistency we
used only the `SPOC` flux and focused on planets with at least three transits, so that we
could make use of the statistical validation tools developed in
[Robnik et al. (2024)](https://arxiv.org/abs/2407.17565). This limited us to consider
$\sim 80$ % of the data, or about 5,000 stars. We report for each planetary candidate the $p$-value of its SNR against the star's own
null distribution in the planet's period bin (the Null Signal Template (NST) test of the
aforementioned paper, with a hierarchical prior pooled over all stars), multiplied by the
star's number of period bins to account for the search over all of them (the global
$p$-value; the local one is listed too, see the Overview).

By default the table shows only the planets that pass every automatic test: no transit better
explained by a false-alarm scenario or by an overlapping deeper planet of the same star, a signal
spread across its transits rather than carried by one, at least three transits, a folded dip shaped
like a transit rather than like a sinusoid at the same period, a SNR significant against the star's
own null distribution, and no known eclipsing binary on the star. The **Show all catalog entries**
button reveals the rest with the failing tests named. The tests are deterministic and recomputed
from the current fits; they are not a human judgement.

Their thresholds are calibrated on an injection and recovery run of the same pipeline: synthetic
planets of 2, 4 and 8 Earth radii are added to real light curves, and each threshold is the one that
recovers the most injected planets at a fixed false alarm rate. Rows whose transit windows contain no
data, or whose transits are covered and show nothing, are not listed here at all.
