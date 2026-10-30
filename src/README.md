# src

Source code for the cryo-jt project, organised by purpose:

- [`data_acquisition`](data_acquisition): code that runs the data
  acquisition (Modbus temperature module, serial pressure transducers, gas
  analyser, real-time plotting). Runs on the lab PC during measurements.
- [`data_analysis`](data_analysis): extraction of the steady-state p-T
  pairs, calculation of Joule-Thomson coefficients with their uncertainties
  (Monte Carlo, uncorrected biases, and method uncertainty), the correction
  and method-uncertainty estimates, the LaTeX tables, and the result figures.
- [`theory_figures`](theory_figures): scripts that regenerate the
  theoretical paper figures (μ_JT-Z, inversion curves, literature
  comparison) from REFPROP, and the measurement-method schematics.

All figure scripts load the shared matplotlib style sheet
[`jced.mplstyle`](jced.mplstyle), which renders text with Helvetica through
LaTeX (packages `helvet` and `newtxsf`), and write into the repository
[`img/`](../img) directory. Each figure is drawn at its printed size (3.3 in
for one column, up to 6.5 in for full width) and is meant to be included
without scaling, so the font sizes in `jced.mplstyle` are the sizes on the
page.

The analysis scripts locate the repository root from their own file location
([`data_analysis/repo_root.py`](data_analysis/repo_root.py)), so they run from
a git clone or from the unpacked Zenodo archive alike; git is not needed.
