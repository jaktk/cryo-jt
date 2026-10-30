# `img/`

Figures used in the paper and the repository documentation, written by the
figure scripts in [`../src/data_analysis/`](../src/data_analysis) and
[`../src/theory_figures/`](../src/theory_figures). All figures use the shared
style sheet [`../src/jced.mplstyle`](../src/jced.mplstyle) (Helvetica through
LaTeX, so the LaTeX packages `helvet` and `newtxsf` are needed) and are drawn
at their printed size.

| Figure | Made by |
| --- | --- |
| `measurements_pures_errorbars.pdf`, `measurements_HeNe_errorbars.pdf`, `measurements_HeN2_errorbars.pdf` | `src/data_analysis/plot_paper_results.py` |
| `theoretical_uncertainty_combined.pdf` | `src/data_analysis/theoretical_jt_uncertainty.py` |
| `theory_vs_experiment_uncertainty.pdf` | `src/data_analysis/compare_theory_vs_experiment.py` |
| `toc.pdf` (table-of-contents graphic, 3.25 in × 1.75 in) | `src/data_analysis/plot_toc_graphic.py` |
| `mu_Z_He_Ne.pdf` | `src/theory_figures/plot_mu_Z_combined.py` (REFPROP) |
| `literature_review_plus_this_work.pdf` | `src/theory_figures/plot_literature_comparison.py` (REFPROP) |
| `JT_theory_small_dp.pdf`, `JT_theory_large_dp.pdf` | `src/theory_figures/plot_jt_method_schematic.py` |
| `Tinv_mix.pdf` | `src/theory_figures/plot_inversion_curves.py` (REFPROP), then finished by hand in Inkscape |
| `PFD.pdf` (simplified flow diagram) | drawing, no script |
| `PID.png`, `photo_*.png` | P&ID drawing and setup photographs, no script |
| `doi-badge-zenodo.svg`, `license-badge-ccby.svg` | README badges |

Regenerating the scripted figures:

```bash
python src/data_analysis/plot_paper_results.py
python src/data_analysis/compare_theory_vs_experiment.py
python src/data_analysis/theoretical_jt_uncertainty.py
python src/data_analysis/plot_toc_graphic.py
python src/theory_figures/plot_jt_method_schematic.py
python src/theory_figures/plot_literature_comparison.py
python src/theory_figures/plot_mu_Z_combined.py
```

`src/theory_figures/plot_inversion_curves.py` is left out on purpose: the
`Tinv_mix.pdf` in this directory was finished by hand in Inkscape (spurious
low-temperature branches removed), and running the script overwrites it with
the unedited plot. Back up the file before rerunning the script.
