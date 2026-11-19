Transit signals found by our search that are **not** in the TESS Objects of Interest (TOI)
catalog. A signal is listed when it

- comes from the blind search of the 2026 runs and has a global null-search p-value below 0.01
  (the null is the star's own noise in the signal's period bin, and the p-value is multiplied by
  the star's number of period bins to account for the search over all of them; see the Overview);
- passes the two automatic vetting flags;
- is not a leftover of a known planet of the same star: its significant transits do not fall on a
  known planet's transit or secondary-eclipse times, and its period is not a known one;
- is not a period duplicate of a stronger candidate.

Every candidate was then inspected by eye and checked against the TOI and community TOI lists,
the NASA Exoplanet Archive and catalogs of eclipsing binaries and variable stars. Signals that
turned out to be known planets or TOIs moved to the catalog, and signals from a neighbouring
binary or variable star were removed.

**Tests.** By default the table shows only the candidates that pass every vetting test: the
pipeline's own tests (spurious transits, per-transit SNR consistency, at least three transits, no
instrumental frequency line), no known eclipsing binary on the star, no match to a single-transit
TOI already listed for the star, no harmonic of a stronger signal of the same star, a folded dip
shaped like a transit rather than like a sinusoid at the same period, and a signal spread across
its transits rather than carried by one. **Show all found candidates** reveals the rest, with the
failed tests listed in the extended columns.

The visual review is not a test: its verdict, confidence and note are columns, and the
`Interesting` flag marks the candidates judged convincing.

**Columns.** *Proposed TOI* continues the host's ExoFOP numbering (TOI 10.01 known, a new
candidate becomes 10.02). *Habitable zone* says whether the candidate's instellation falls in the
conservative or the optimistic habitable zone of
[Kopparapu et al. (2014)](https://doi.org/10.1088/2041-8205/787/2/L29) (see the Overview).

**Interactive plot** in the toolbar shows the population in any pair of catalog quantities.
Convincing candidates are drawn as orange stars, and clicking a point opens its figures.
