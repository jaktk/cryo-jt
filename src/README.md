# src

Source code for the cryo-jt project, organised by purpose:

- [`data_acquisition`](data_acquisition) — code that runs the data
  acquisition (Modbus temperature module, serial pressure transducers, gas
  analyser, real-time plotting). Runs on the lab PC during measurements.
- [`data_analysis`](data_analysis) — data wrangling, calculation of
  Joule-Thomson coefficients with conventional and Monte Carlo
  uncertainties, and the publication-style result figures. Shares the
  matplotlib style sheet
  [`jced.mplstyle`](jced.mplstyle) used by every figure script.
- [`theory_figures`](theory_figures) — scripts that regenerate the
  theoretical paper figures (μ_JT–Z, inversion curves, literature
  comparison) from REFPROP, in the same style.

Figure scripts write into the repository [`img/`](../img) directory.
