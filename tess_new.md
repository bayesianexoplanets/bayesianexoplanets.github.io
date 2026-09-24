New planet candidates found by our pipeline that are **not** in the TESS Objects of
Interest (TOI) catalog: every blind-search candidate of the 2026 searches with a
global null-simulation p-value below 0.01 (the null is the star's own noise in the
candidate's period bin, and the p-value is multiplied by the star's number of period
bins to account for the search over all of them, see the Overview), passing the two automatic
vetting flags, not a leftover of a known planet of the same star (its significant transits do not
fall on a known planet's transit or secondary-eclipse times, and its period is not a known one), and
not a period duplicate of a stronger candidate. Every candidate was inspected by eye and checked
against the TOI and community TOI lists, the NASA Exoplanet Archive and catalogs of eclipsing
binaries and variable stars: signals that turned out to be known planets or TOIs moved to the
catalog, and signals from a neighbouring binary or variable star were removed. By default the table shows only the candidates that pass every vetting test: the
pipeline's own tests (spurious transits, per-transit SNR consistency, at least three
transits, no instrumental frequency line), no known eclipsing binary on the star, no
match to a single-transit TOI already listed for the star, no harmonic of a
stronger signal of the same star, a folded dip shaped
like a transit rather than like a sinusoid at the same period, and a signal spread across its
transits rather than carried by one. The visual review is not a test: its verdict, confidence and
note are columns, and the `Interesting` flag marks the candidates judged convincing. The
**Show all found candidates** button reveals the rest, with the failed tests listed in
the extended columns.

Use the **Interactive plot** button in the toolbar to explore the population in any
pair of catalog quantities. Convincing candidates are drawn as orange stars, and
clicking any point opens its diagnostic plots.
