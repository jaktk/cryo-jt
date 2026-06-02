from __future__ import annotations

import argparse
import logging
import os
import sys
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from numpy.polynomial.chebyshev import Chebyshev

from FluidProps import FluidProps
from get_git_root import get_git_root
from theoretical_jt_uncertainty import U_P, combined_temperature_uncertainty

logger = logging.getLogger(__name__)

DEFAULT_N_SAMPLES = 2000
DEFAULT_SEED = 20210101
DEFAULT_U_X_K2 = 1e-3  # gas analyser expanded uncertainty, 0.1 mol-% (k=1.96)
COVERAGE_FACTOR = 1.96  # 95 %

@dataclass(frozen=True)
class IsenthalpInputs:
    filename: str
    fluid: str | tuple[str, str]
    is_mixture: bool
    polynomial_degree: int
    p_in_MPa: float
    T_in_K: float
    p_measured_MPa: np.ndarray
    x1_mean: float | None # first component mole fraction; None for pure fluids


def monte_carlo_bound(inputs: IsenthalpInputs,
                      n_samples: int = DEFAULT_N_SAMPLES,
                      u_x_k2: float = DEFAULT_U_X_K2,
                      rng: np.random.Generator | None = None,) -> pd.DataFrame:

    if rng is None:
        rng = np.random.default_rng(DEFAULT_SEED)

    fp = _make_fluid_props(inputs)

    # Build the idealised isenthalp: place the measured pressures on the
    # EOS h=const. line through the mean inlet condition.
    h_in = fp.get_hmass(p=inputs.p_in_MPa, T=inputs.T_in_K)
    p_meas = inputs.p_measured_MPa
    T_eos = np.array([fp.get_T(p=p, h=h_in) for p in p_meas])
    jt_eos = np.array([fp.get_JT_coefficient(p=p, T=T) for p, T in zip(p_meas, T_eos)])

    # Sensor-floor uncertainties at the idealised points (k=1.96).
    u_T = combined_temperature_uncertainty(T_eos)  # K, k=1.96 expanded
    u_p = np.full_like(p_meas, U_P, dtype=float)  # MPa, k=1.96 expanded
    # 1-sigma half-widths for the uniform distribution
    # (uniform on [-U, +U] has std = U/sqrt(3); we draw on [-U, +U] to match
    # the paper's "random errors within the bounds of Table 3" convention).
    n_points = p_meas.size
    deg = max(1, min(inputs.polynomial_degree, n_points - 1))

    jt_samples = np.empty((n_samples, n_points), dtype=float)
    jt_samples.fill(np.nan)
    n_failed = 0

    for k in range(n_samples):
        p_pert = p_meas + rng.uniform(-u_p, u_p)
        T_pert = T_eos + rng.uniform(-u_T, u_T)
        try:
            fit = Chebyshev.fit(p_pert, T_pert, deg)
            deriv = fit.deriv(m=1)
            jt_samples[k, :] = deriv(p_meas)
        except (np.linalg.LinAlgError, ValueError):
            n_failed += 1
            continue

    if n_failed:
        logger.warning(
            "%s: %d of %d Chebyshev fits failed and were dropped",
            inputs.filename, n_failed, n_samples,
        )

    valid = ~np.isnan(jt_samples).any(axis=1)
    if not valid.any():
        raise RuntimeError(
            f"{inputs.filename}: every Monte Carlo realisation failed to fit; "
            "check the polynomial degree and inlet conditions"
        )

    jt_samples = jt_samples[valid]
    jt_mc_mean = jt_samples.mean(axis=0)
    jt_mc_std = jt_samples.std(axis=0, ddof=1)
    jt_mc_unc = COVERAGE_FACTOR * jt_mc_std

    with np.errstate(divide="ignore", invalid="ignore"):
        jt_mc_rel_unc_perc = 100.0 * jt_mc_unc / np.abs(jt_eos)

    df = pd.DataFrame({
        "p/MPa": p_meas,
        "T_eos/K": T_eos,
        "jt_eos/(K/MPa)": jt_eos,
        "jt_mc_mean/(K/MPa)": jt_mc_mean,
        "jt_mc_std/(K/MPa)": jt_mc_std,
        "jt_mc_unc_k2/(K/MPa)": jt_mc_unc,
        "jt_mc_rel_unc_perc": jt_mc_rel_unc_perc,
    })

    if inputs.is_mixture:
        df = _add_composition_contribution(
            df, inputs, fp, p_meas, T_eos, jt_eos, deg, u_p, u_T,
            n_samples, u_x_k2, rng,
        )
    return df


def _make_fluid_props(inputs: IsenthalpInputs) -> FluidProps:
    fp = FluidProps(inputs.fluid)
    if inputs.is_mixture:
        if inputs.x1_mean is None:
            raise ValueError(
                f"{inputs.filename}: mixture has no composition recorded"
            )
        fp.set_composition_from_1st_fraction(inputs.x1_mean)
    return fp


def _add_composition_contribution(base: pd.DataFrame,
                                  inputs: IsenthalpInputs,
                                  fp: FluidProps,
                                  p_meas: np.ndarray,
                                  T_eos: np.ndarray,
                                  jt_eos: np.ndarray,
                                  deg: int,
                                  u_p: np.ndarray,
                                  u_T: np.ndarray,
                                  n_samples: int,
                                  u_x_k2: float,
                                  rng: np.random.Generator) -> pd.DataFrame:
    n_points = p_meas.size
    jt_samples = np.full((n_samples, n_points), np.nan)
    n_failed = 0

    x1 = float(inputs.x1_mean)  # type: ignore[arg-type]
    # Clip the uniform sampling range so x1 stays in (0, 1)
    low = max(-u_x_k2, -x1)
    high = min(u_x_k2, 1 - x1)
    # EOS lookups near a phase boundary occasionally return non-physical
    # values when composition is perturbed; reject realisations where any
    # per-point JT differs from the nominal-composition EOS by more than 5x
    max_rel_eos_deviation = 5.0

    # Composition-only contribution: perturb the mixture composition alone while holding p and T on the EOS isenthalp
    for k in range(n_samples):
        x_pert = x1 + rng.uniform(low, high, size=n_points)
        try:
            jt_calc = np.empty(n_points)
            for i, (p, T, x) in enumerate(zip(p_meas, T_eos, x_pert)):
                fp.set_composition_from_1st_fraction(x)
                jt_calc[i] = fp.get_JT_coefficient(p=p, T=T)
            if not np.all(np.isfinite(jt_calc)):
                n_failed += 1
                continue
            if np.any(np.abs(jt_calc - jt_eos) > max_rel_eos_deviation * np.abs(jt_eos)):
                n_failed += 1
                continue
            jt_samples[k, :] = jt_calc
        except (np.linalg.LinAlgError, ValueError):
            n_failed += 1

    # Restore the nominal composition for the FluidProps caller
    fp.set_composition_from_1st_fraction(x1)

    if n_failed:
        logger.warning(
            "%s: %d of %d composition-MC fits failed and were dropped",
            inputs.filename, n_failed, n_samples,
        )
    valid = ~np.isnan(jt_samples).any(axis=1)
    if not valid.any():
        logger.warning(
            "%s: composition contribution could not be evaluated", inputs.filename
        )
        base["jt_mc_unc_cmp_k2_perc"] = np.nan
        return base

    jt_samples = jt_samples[valid]
    jt_mc_unc_cmp = COVERAGE_FACTOR * jt_samples.std(axis=0, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        base["jt_mc_unc_cmp_k2_perc"] = 100.0 * jt_mc_unc_cmp / np.abs(jt_eos)
    return base


def load_isenthalp(csv_path: str, fluid_name: str, polynomial_degree: int) -> IsenthalpInputs:
    df = pd.read_csv(csv_path)
    fluid, is_mixture = _parse_fluid(fluid_name)

    p_meas = df["PT102/MPa"].to_numpy()
    p_in = float(df["PT101/MPa"].mean())
    T_in = float(df["TT101/K"].mean())

    x1_mean: float | None = None
    if is_mixture and "x1" in df.columns:
        x1 = df["x1"].dropna()
        x1_mean = float(x1.mean()) if not x1.empty else None

    return IsenthalpInputs(
        filename=os.path.basename(csv_path),
        fluid=fluid,
        is_mixture=is_mixture,
        polynomial_degree=int(polynomial_degree),
        p_in_MPa=p_in,
        T_in_K=T_in,
        p_measured_MPa=p_meas,
        x1_mean=x1_mean,
    )


def _parse_fluid(name: str) -> tuple[str | tuple[str, str], bool]:
    parts = name.split("-")
    if len(parts) == 1:
        return parts[0], False
    if len(parts) == 2:
        return (parts[0], parts[1]), True
    raise ValueError(f"unexpected fluid name {name!r}")


def _load_summary(repo_root: str) -> pd.DataFrame:
    path = os.path.join(repo_root, "data", "derived_data", "jt_coeffs",
                        "jt_coefficients_summary.csv")
    return pd.read_csv(path)


def run_all(repo_root: str,
            output_dir: str,
            n_samples: int,
            u_x_k2: float,
            rng: np.random.Generator,
            only: Iterable[str] | None = None) -> pd.DataFrame:
    
    """Run the Monte Carlo for every isenthalp in the summary file."""
    os.makedirs(output_dir, exist_ok=True)
    summary = _load_summary(repo_root)
    pt_dir = os.path.join(repo_root, "data", "derived_data", "p_T_pairs")
    only_set = {os.path.splitext(s)[0] for s in only} if only else None

    rows: list[dict] = []
    for _, meta in summary.iterrows():
        filename: str = meta["filename"]
        stem = os.path.splitext(filename)[0]
        if only_set is not None and stem not in only_set:
            continue
        csv_path = os.path.join(pt_dir, filename)
        if not os.path.exists(csv_path):
            logger.warning("missing p_T_pairs CSV for %s; skipping", filename)
            continue

        inputs = load_isenthalp(
            csv_path,
            fluid_name=meta["fluid"],
            polynomial_degree=int(meta["polynomial_degree"]),
        )
        logger.info(
            "Monte Carlo: %s (n=%d, deg=%d, mixture=%s)",
            stem, inputs.p_measured_MPa.size, inputs.polynomial_degree,
            inputs.is_mixture,
        )
        df = monte_carlo_bound(inputs, n_samples=n_samples, u_x_k2=u_x_k2, rng=rng)
        df.to_csv(os.path.join(output_dir, f"{stem}_MC_uncertainty.csv"), index=False)

        rel = df["jt_mc_rel_unc_perc"]
        row = {
            "filename": filename,
            "fluid": meta["fluid"],
            "is_mixture": bool(meta["is_mixture"]),
            "x1": inputs.x1_mean,
            "n_points": int(inputs.p_measured_MPa.size),
            "polynomial_degree": int(inputs.polynomial_degree),
            "mean_jt_mc_rel_unc_perc": float(rel.mean()),
            "max_jt_mc_rel_unc_perc": float(rel.max()),
        }
        if "jt_mc_unc_cmp_k2_perc" in df.columns:
            cmp = df["jt_mc_unc_cmp_k2_perc"]
            row["mean_jt_mc_cmp_rel_unc_perc"] = float(cmp.mean())
            row["max_jt_mc_cmp_rel_unc_perc"] = float(cmp.max())
        rows.append(row)

    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(output_dir, "mc_uncertainty_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info("wrote summary to %s (%d isenthalps)", summary_path, len(summary_df))
    return summary_df


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Monte Carlo bound on the achievable mu_JT uncertainty.",
    )
    parser.add_argument(
        "--n-samples", type=int, default=DEFAULT_N_SAMPLES,
        help="number of Monte Carlo realisations per isenthalp",
    )
    parser.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help="seed for the NumPy random generator",
    )
    parser.add_argument(
        "--u-x-k2", type=float, default=DEFAULT_U_X_K2,
        help="expanded (k=1.96) composition uncertainty, mole fraction",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help=("directory for the per-isenthalp CSVs and the summary; "
              "defaults to data/derived_data/jt_coeffs/"),
    )
    parser.add_argument(
        "--isenthalp", action="append", default=None,
        metavar="STEM",
        help=("isenthalp stem(s) to run, e.g. Helium-Neon_65K_5MPa; "
              "may be passed multiple times. Default: all"),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-q", "--quiet", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    level = logging.DEBUG if args.verbose else logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s %(name)s: %(message)s")

    repo_root = get_git_root(os.getcwd())
    output_dir = args.output_dir or os.path.join(
        repo_root, "data", "derived_data", "jt_coeffs",
    )
    rng = np.random.default_rng(args.seed)
    run_all(
        repo_root=repo_root,
        output_dir=output_dir,
        n_samples=args.n_samples,
        u_x_k2=args.u_x_k2,
        rng=rng,
        only=args.isenthalp,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
