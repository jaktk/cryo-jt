# `data/derived_data/jt_coeffs/`

Derived Joule-Thomson coefficients with uncertainties, one isenthalp per
file, plus a summary table over all isenthalps and the same content as JSON
for programmatic use. The Monte Carlo bound on the achievable uncertainty
(`*_MC_uncertainty.csv`) is written alongside as a separate set of files.

These products are computed from the steady-state $p$-$T$ pairs in
[`../p_T_pairs/`](../p_T_pairs) by
[`../../../src/data_analysis/calculate_jt_coefficients.py`](../../../src/data_analysis/calculate_jt_coefficients.py),
with the Monte Carlo bound added by
[`../../../src/data_analysis/monte_carlo_uncertainty.py`](../../../src/data_analysis/monte_carlo_uncertainty.py).
The full column dictionary is in
[`../../metadata/codebook.csv`](../../metadata/codebook.csv).

## Files

| File | Content |
| --- | --- |
| `Fluid_TK_PMPa_JT.csv` | One isenthalp: per-point pressure, temperature, derived $\mu_{\rm JT}$ with conventional and Monte Carlo expanded uncertainties (k=2), and the EOS reference at each point. |
| `Fluid_TK_PMPa_MC_uncertainty.csv` | One isenthalp: the Monte Carlo bound on the achievable $\mu_{\rm JT}$ uncertainty given the sensor floor (and, for mixtures, the composition uncertainty), evaluated on the EOS-idealised isenthalp through the mean inlet condition. |
| `jt_coefficients_summary.csv` | One row per isenthalp: fluid, composition, fit metadata, mean conventional uncertainty, and mean deviation from the EOS. |
| `jt_coefficients_detailed.json` | The same results in a nested JSON form including the per-point arrays, for programmatic use. |
| `mc_uncertainty_summary.csv` | One row per isenthalp: mean and maximum Monte Carlo relative expanded uncertainty (and, for mixtures, the marginal composition contribution). |

## Filename convention

`<Fluid>_<T>K_<P>MPa_JT.csv` and `<Fluid>_<T>K_<P>MPa_MC_uncertainty.csv`,
where:

| Token | Meaning |
| --- | --- |
| `Fluid` | `Nitrogen`, `Argon`, `Helium`, `Helium-Neon`, `Nitrogen-Helium` (for mixtures the order matches the `x1`/`x2` columns) |
| `T` | nominal inlet temperature, rounded to K |
| `P` | nominal inlet pressure, rounded to MPa |

The same stem is used for the matching p-T input in
[`../p_T_pairs/`](../p_T_pairs).

## Uncertainty conventions

Both the conventional (GUM) and the Monte Carlo uncertainty are reported.
The conventional one (Eq.~9 in the paper) is propagated in quadrature from
the sensor-chain uncertainties of Table~3 and gives `mean_unc_perc` in the
summary file. The Monte Carlo one (1000 samples around each measured (p, T)
within its sensor uncertainty, refitted with the same Chebyshev polynomial
degree) is written to `JT_UNC/(K/MPa)` in the per-isenthalp files at 1σ;
the publication figures plot it at k=2 (coverage factor 2, 95.45 % level
of confidence).
