# Theory Figures

Scripts that generate the *theoretical* (data-independent) figures of the
paper directly from REFPROP / CoolProp through the shared
[`FluidProps`](../data_analysis/FluidProps.py) wrapper, in the
same matplotlib style sheet ([`jced.mplstyle`](../jced.mplstyle))
as the result figures.

All scripts write into the repository [`img/`](../../img) directory.

## Scripts

| Script | Output(s) in `img/` |
| --- | --- |
| [`plot_mu_Z_combined.py`](plot_mu_Z_combined.py) | `mu_Z_He_Ne.pdf` |
| [`plot_inversion_curves.py`](plot_inversion_curves.py) | `Tinv_mix.pdf` |
| [`plot_literature_comparison.py`](plot_literature_comparison.py) | `literature_review_plus_this_work.pdf` |

### `plot_mu_Z_combined.py`

Two-panel μ_JT and Z (compressibility factor) versus temperature for pure
helium-4 and pure neon, with pressure shown by a viridis colour bar.

### `plot_inversion_curves.py`

⁴He-Ne inversion curves in p-T space at composition
x_He ∈ {0.0, 0.1, …, 1.0} (viridis colour bar), with the p-T region
targeted by the measurements of this work shaded for reference. Writes
`Tinv_mix.pdf`. Each inversion pressure is obtained by a Brent root-find on
μ_JT(p, T) = 0 in [0.05, 45] MPa.

### `plot_literature_comparison.py`

Relative deviation of historical Joule-Thomson measurements (and the
present work) from the corresponding accurate Helmholtz-energy EOS, plotted
versus temperature. Each point is coloured by author and shaped by fluid; a
±5 % band is shaded. The script reads
[`data/literature/literature_measurements.csv`](../../data/literature/literature_measurements.csv)
and evaluates the EOS μ_JT via REFPROP.

## Requirements

These scripts require [REFPROP](https://www.nist.gov/srd/refprop) via the
[`FluidProps`](../data_analysis/FluidProps.py) wrapper, plus
`numpy`, `pandas`, `scipy`, and `matplotlib`. The shared style sheet uses
`text.usetex` and therefore needs a LaTeX toolchain (`latex` + `dvipng`); to
regenerate without LaTeX, set `text.usetex` to `False` in
[`jced.mplstyle`](../jced.mplstyle).

Regenerate the full theory-figure set with:

```bash
python src/theory_figures/plot_mu_Z_combined.py
python src/theory_figures/plot_inversion_curves.py
python src/theory_figures/plot_literature_comparison.py
```
