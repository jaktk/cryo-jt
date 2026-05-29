# Literature data

[`literature_measurements.csv`](literature_measurements.csv) — compilation of
historical Joule-Thomson coefficient measurements from the literature
(Bier *et al.* 1974, King & Potter 1962, Roebuck & Osterberg 1933/1934/1935,
Roebuck 1926, Smith 1970) together with the present-work measurements
(Tkaczuk). Used by
[`src/theory_figures/plot_literature_comparison.py`](../../src/theory_figures/plot_literature_comparison.py)
to build the `literature_review_plus_this_work.pdf` figure.

Columns: `authors, year, fluid, p/psia, p/atm, p/MPa, T/degC, T/F, T/K, x1,
x2, uJT_meas/(K/MPa), uJT_meas/(F/psi), uJT_meas/(degC/atm), rel_err`.

Fluids covered: propylene, n-butane, nitrogen, argon, air, helium,
helium-neon, helium-nitrogen.
