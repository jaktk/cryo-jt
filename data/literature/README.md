# Literature data

[`literature_measurements.csv`](literature_measurements.csv): compilation of
historical Joule-Thomson coefficient measurements from the literature
(Bier *et al.* 1974, King & Potter 1962, Roebuck & Osterberg 1933/1934/1935,
Roebuck 1926, Smith 1970). Used by
[`src/theory_figures/plot_literature_comparison.py`](../../src/theory_figures/plot_literature_comparison.py)
to build the `literature_review_plus_this_work.pdf` figure. The points of this
work in that figure are read from the pipeline output
(`data/derived_data/jt_coeffs/*_JT.csv`, points with `used` true), so they
always match the result tables.

Up to version 1.3.1 the file also held 74 rows labelled `Tkaczuk,2021` (22 N2,
6 Ar, 6 He, 28 He-Ne, 12 He-N2) from an earlier analysis of the same
measurements. They did not match the current pipeline output and were removed
in version 1.3.2.

Columns: `authors, year, fluid, p/psia, p/atm, p/MPa, T/degC, T/F, T/K, x1,
x2, uJT_meas/(K/MPa), uJT_meas/(F/psi), uJT_meas/(degC/atm), rel_err`.

Fluids covered: propylene, n-butane, nitrogen, argon, air, helium.
