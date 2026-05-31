# `img/`

Figures used in the paper and the repository documentation, written by the
figure scripts in [`../src/data_analysis/`](../src/data_analysis) and
[`../src/theory_figures/`](../src/theory_figures). The Piping and
Instrumentation Diagram (`PID.png`) and the setup photographs are
hand-prepared and committed directly.

Each script and the figure it produces are listed in
[`../src/data_analysis/README.md`](../src/data_analysis/README.md).
Regenerating the full figure set:

```bash
python src/data_analysis/plot_paper_results.py
python src/data_analysis/compare_theory_vs_experiment.py
python src/data_analysis/theoretical_jt_uncertainty.py
python src/data_analysis/plot_toc_graphic.py
python src/theory_figures/plot_inversion_curves.py
python src/theory_figures/plot_literature_comparison.py
python src/theory_figures/plot_mu_Z_combined.py
```

The schematic method figures (`JT_theory_small_dp.pdf`,
`JT_theory_large_dp.pdf`), the simplified flow diagram (`PFD.pdf`), the P&ID
(`PID.png`), and the setup photographs are hand-prepared and committed
directly; they have no generating script. `nonhTransform.pdf` is an auxiliary
sketch not used in the paper.
