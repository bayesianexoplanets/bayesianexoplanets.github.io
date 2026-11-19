Eclipsing binaries on TESS stars, collected so that a transit search does not rediscover them.
The list combines three sources (**Listed as**):

- *TESS-EB catalog*: the TESS eclipsing-binary catalog of
  [Prša et al. (2022)](https://doi.org/10.3847/1538-4365/ac324a), 4,476 entries;
- *TOI disposition*: TESS Objects of Interest that the follow-up working group dispositioned as
  binaries, 639;
- *false-positive analysis* and *review*: planet candidates that we reclassified. For 647 the
  binary scenarios of a Bayesian comparison of the folded light curve exceed 50 % probability
  (**P(binary)**); 10 show a secondary eclipse or alternating eclipse depths in the review.

**Checked.** For every binary in our archive the period, epoch and eclipse width were compared
with the light curves, sector by sector (see *Known planets and binaries* in the Overview). 5,135
show the eclipse (*clear* 4,671, *weak* 464) and carry the values read off the data.

A *second review* went back to 325 binaries with the gentler detrending of the figures below: 245
where the first review had seen no eclipse because its detrending had absorbed it, and 80 whose
fitted curve missed the eclipse. It found the eclipse in 315, moved the epoch or the period of 218
(half of the epoch shifts exceed 2.6 hours, and a tenth exceed two days), and widened most
eclipses: the listed widths were too narrow by a median factor 1.5. Five signals are periodic but
not eclipses. Hidden test rows with a deliberately wrong epoch or period measure how reliable this
is: after two passes 34 of 36 were right, and both misses were periods off by a cycle count
across a gap of years. So about one in twenty of these 315 ephemerides may still be wrong, most
likely where the period was changed.

203 binaries show no eclipse in either review and keep the catalog values, which is not evidence
against the binary. 368 rest on the false-positive analysis alone, and 61 have not been checked.

**Columns.** *Eclipse width* is the full width of the primary eclipse, from first to last contact
(for V-shaped eclipses the full base). *Secondary phase* is the position of the secondary eclipse
as a fraction of the period after the primary, from the catalog or measured in the review; 4,354
binaries have one. *SNR* exists for the 718 binaries that our transit search had fitted as planet
candidates.

**Figures.** 5,750 binaries have a folded light curve, from every TESS sector of the star in our
archive: the primary and the secondary eclipse on top and the whole orbit below; contact binaries
show the whole orbit only. The light curve is detrended in two steps, so that the detrending
follows the star's variability and leaves the eclipses alone. **Curve** says what the dashed line
is:

- *posterior* (646): for reclassified candidates, the binary model of the false-positive analysis
  itself, the median of its posterior samples with a 68 % band. Only a baseline and a time shift
  are fitted, and the fold uses the two-minute sectors that analysis was fitted to. Where the
  curve misses the data, that is a statement about the model. 209 of these are pairs of similar
  stars at twice the listed period;
- *physical fit* (2,864): an eclipse model of two limb-darkened stars fitted to the fold;
- *two-Gaussian fit* (2,016): two Gaussian eclipses and a sinusoid fitted to the fold, for contact
  binaries and where it describes the fold better. The parameters of both fits are degenerate
  and are not measurements, so we do not list them;
- *none* (203): no review found an eclipse at the listed place; the figure shows the data only;
- *no curve* (21): the eclipse is there, but a fitted curve does not describe it, mostly a narrow
  eclipse on top of a larger wave; the figure shows the data only.

A reclassified candidate is a statement about probabilities, not a confirmation: 647 of these
rows are binaries only by the model comparison, and a planet whose catalogued period is twice
the true one can look like a binary with equal eclipses.
