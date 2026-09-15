New planet candidates found by our pipeline that are **not** in the TESS Objects of
Interest (TOI) catalog: every blind-search candidate of the 2026 rerun with a
period-local null-simulation p-value below $10^{-4}$ (the null is the star's own
noise in the candidate's period bin, see the Overview), passing the two automatic
vetting flags, not a harmonic of a known planet of the same star, and not a period
duplicate of a stronger candidate. By default the table shows only the candidates that pass every vetting test: the
pipeline's own tests (spurious transits, per-transit SNR consistency, at least three
transits, no instrumental frequency line), no known eclipsing binary on the star, no
match to a single-transit TOI already listed for the star, no harmonic of a
stronger signal of the same star, and no integer or half-integer period ratio with a
known planet of the star. The visual review is not a test: its verdict, confidence and
note are columns, and the `Interesting` flag marks the candidates judged convincing. The
**Show all found candidates** button reveals the rest, with the failed tests listed in
the extended columns.

Use the **Interactive plot** button in the toolbar to explore the population in any
pair of catalog quantities. Convincing candidates are drawn as orange stars, and
clicking any point opens its diagnostic plots.
