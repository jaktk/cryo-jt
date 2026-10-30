"""
Weighted Chebyshev fit of T(p) along an isenthalp and selection of its degree.

Shared by calculate_jt_coefficients.py, monte_carlo_uncertainty.py,
extract_fit_coefficients.py, and validate_degree_criterion.py, so that every
script fits the same way. No REFPROP needed.

Degree rule: the lowest degree n (1 to MAX_DEGREE, at most N_s - 2) whose
weighted residuals agree with the measurement uncertainty, i.e. whose reduced
chi-square is below the CHI2_LEVEL quantile of chi2_nu / nu with nu = N_s - n - 1
degrees of freedom. If no degree up to the cap N_s - 2 passes, the degree with
the lowest reduced chi-square is used and the fit is marked as not passing.
The cap keeps at least one residual degree of freedom, so the polynomial never
interpolates. The only exception is a three-point
isenthalp, which is interpolated with n = 2: a straight line through three
unevenly spaced points biases the slope at the middle point.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.chebyshev import Chebyshev
from scipy import stats

MAX_DEGREE = 6
CHI2_LEVEL = 0.95


@dataclass
class FitResult:
    polynomial: Chebyshev
    degree: int
    dof: int            # residual degrees of freedom, N_s - n - 1
    chi2_red: float     # reduced chi-square chi2 / dof (NaN when dof = 0)
    chi2_crit: float    # CHI2_LEVEL quantile of chi2_dof / dof (NaN when dof = 0)
    passed: bool | None # chi2_red <= chi2_crit; None for an interpolation
    s_T: float          # RMS of the temperature residuals, K


def effective_t_uncertainty(pressures, temperatures, u_p, u_T):
    """
    Standard uncertainty of each point expressed in temperature: the pressure
    uncertainty is folded in through the local slope dT/dp, estimated by finite
    differences, u_eff = sqrt(u_T^2 + (u_p * dT/dp)^2). Inputs and output share
    one coverage factor (standard uncertainties for the chi-square test). """
    order = np.argsort(pressures)
    slope = np.empty_like(temperatures, dtype=float)
    slope[order] = np.gradient(temperatures[order], pressures[order])
    return np.sqrt(u_T**2 + (u_p * slope)**2)


def chebyshev_fit_weights(pressures, temperatures, u_p, u_T):
    """
    Weights for Chebyshev.fit. NumPy applies w to the unsquared residual, so
    w = 1 / u_eff weights the squared residuals by 1 / u_eff^2 (inverse-variance
    weighting). """
    return 1.0 / effective_t_uncertainty(pressures, temperatures, u_p, u_T)


def fit_chebyshev(pressures, temperatures, degree, weights=None):
    """ Least-squares Chebyshev fit, returned in the raw-pressure basis
    T/K = sum_k a_k C_k(p/MPa) (SI eq S5). """
    fit = Chebyshev.fit(pressures, temperatures, degree, w=weights)
    return Chebyshev(fit.convert().coef)


def fit_statistics(pressures, temperatures, polynomial, u_eff):
    """ Reduced chi-square, its critical value, and the residual RMS of a fit. """
    residuals = temperatures - polynomial(pressures)
    dof = len(pressures) - polynomial.degree() - 1
    if dof > 0:
        chi2_red = np.sum((residuals / u_eff)**2) / dof
        chi2_crit = stats.chi2.ppf(CHI2_LEVEL, dof) / dof
    else:
        chi2_red = chi2_crit = np.nan
    return dof, chi2_red, chi2_crit, np.sqrt(np.mean(residuals**2))


def select_degree(pressures, temperatures, u_p, u_T, max_degree=MAX_DEGREE):
    """
    Fit an isenthalp with the degree rule of the module docstring.
    u_p and u_T are standard (k = 1) uncertainties of each point.
    """
    n_points = len(pressures)
    if n_points < 3:
        raise ValueError(f"need at least three points, got {n_points}")
    u_eff = effective_t_uncertainty(pressures, temperatures, u_p, u_T)
    weights = 1.0 / u_eff

    if n_points == 3:
        poly = fit_chebyshev(pressures, temperatures, 2, weights)
        dof, chi2_red, chi2_crit, s_t = fit_statistics(pressures, temperatures, poly, u_eff)
        return FitResult(poly, 2, dof, chi2_red, chi2_crit, None, s_t)

    cap = min(max_degree, n_points - 2)
    tried = []
    for degree in range(1, cap + 1):
        poly = fit_chebyshev(pressures, temperatures, degree, weights)
        dof, chi2_red, chi2_crit, s_t = fit_statistics(pressures, temperatures, poly, u_eff)
        if chi2_red <= chi2_crit:
            return FitResult(poly, degree, dof, chi2_red, chi2_crit, True, s_t)
        tried.append(FitResult(poly, degree, dof, chi2_red, chi2_crit, False, s_t))
    # No degree passes: keep the one with the lowest reduced chi-square
    return min(tried, key=lambda f: f.chi2_red)
