import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from numpy.polynomial.chebyshev import Chebyshev

STYLE = os.path.join(os.path.dirname(__file__), '..', 'jced.mplstyle')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'img')

U_P = 0.0001 * 13.7  # MPa => 1.37e-3 MPa # Pressure sensor uncertainty: Mensor CPT 6100, 0.01% FS, p_max = 13.7 MPa

def combined_temperature_uncertainty(T):
    """
    Combined temperature measurement uncertainty in K.
    Sum of Cernox sensor, CABTR acquisition module, and calibration-polynomial
    fit, using the calibration of the outlet thermometer TT102 (sensor X93303),
    which sets the fitted p-T pairs and so matches the data-reduction pipeline.
    Returns uncertainty for coverage factor k=1.96 (95%).

    Temperature in K (valid 20-325 K)
    """
    T = np.asarray(T, dtype=float)
    u_cernox = np.maximum(-5e-05 * T**2 + 0.1281 * T + 6.4603, 0.5) # Cernox sensor uncertainty [mK]
    u_cabtr = np.maximum(2e-05 * T**2 + 0.0695 * T - 0.1001, 0.5) # CABTR acquisition module uncertainty [mK]

    # Calibration-polynomial fit uncertainty [mK], sensor X93303, piecewise over
    # the two calibration bands: 14.1-80 K (N=31, n=7, DTrms=1.60 mK) and
    # 80-325 K (N=32, n=8, DTrms=5.84 mK); expanded to k=2.
    def _poly(t):
        N, n, DTrms = (31, 7, 1.60) if t < 80.0 else (32, 8, 5.84)
        return 2 * (N / (N - n) * DTrms ** 2) ** 0.5
    u_poly = np.vectorize(_poly)(T)
    # The three contributions above are expanded at k = 2; rescale the combined
    # chain to k = 1.96 (95 % level of confidence) used throughout the analysis.
    return (u_cernox + u_cabtr + u_poly) * 1e-3 * (1.96 / 2.0)  # combined [K]


def conventional_relative_uncertainty(mu_jt, delta_T, delta_p, T_in, T_out, k=1.96):
    """
    Expanded relative standard uncertainty of the JT coefficient
    from conventional error propagation (Eq. 5.10 in thesis).

    U_r(mu_JT) / |mu_JT| = k * sqrt([U(T_in)^2 + U(T_out)^2] / (T_in - T_out)^2
                         + [U(p_in)^2 + U(p_out)^2] / (p_in - p_out)^2)
    """
    U_Tin = combined_temperature_uncertainty(T_in)
    U_Tout = combined_temperature_uncertainty(T_out)
    U_pin = U_P
    U_pout = U_P

    term_T = (U_Tin**2 + U_Tout**2) / delta_T**2
    term_p = (U_pin**2 + U_pout**2) / delta_p**2
    return k * np.sqrt(term_T + term_p)


def monte_carlo_uncertainty(slope, p_in, p_out, n_points, T_in=80.0, curvature=0.0, n_iter=1500, k=1.96):
    """
    Monte Carlo uncertainty estimation for the differentiated JT coefficient, including the polynomial fitting step.

    The synthetic isenthalp is parameterized with p_in as the high-pressure inlet (where T = T_in) and p_out as the low-pressure outlet.
    slope = dT/dp > 0 means the fluid cools on expansion (positive mu_JT). The line is:

        T(p) = T_in + slope*(p - p_in) + curvature*(p - p_in)^2

    so T_out = T_in + slope*(p_out - p_in) < T_in for slope > 0 since p_out < p_in.

    Returns
    -------
    mean_rel_err : float
        Mean relative expanded uncertainty (k=1.96) across interior points.
    point_rel_errs : array
        Relative expanded uncertainty at each interior point.
    """
    p_meas = np.linspace(p_in, p_out, n_points) # measurement points equally spaced in pressure (descending)
    T_meas = T_in + slope * (p_meas - p_in) + curvature * (p_meas - p_in)**2 # synthetic isenthalp: T(p) with T_in at the high-pressure end
    T_meas = np.maximum(T_meas, 20.0) # guard against unphysical temperatures for very steep slopes
    
    mu_true = slope + 2.0 * curvature * (p_meas - p_in) # true JT coefficient at each point (analytical derivative dT/dp)

    U_T = np.array([combined_temperature_uncertainty(T) for T in T_meas]) # temperature uncertainty at each point (1-sigma = U_T/k)
    sigma_T = U_T / k
    sigma_p = U_P / k

    # Choose polynomial degree (set by the isenthalp shape)
    # Linear isenthalps need degree 2, curved isenthalps need degree 3.
    degree = min(3, n_points - 1) if abs(curvature) > 1e-10 else min(2, n_points - 1)
    degree = max(1, degree)

    # Monte Carlo iterations
    mu_samples = np.zeros((n_iter, n_points))
    for i in range(n_iter):
        p_pert = p_meas + np.random.normal(0, sigma_p, n_points)
        T_pert = T_meas + np.random.normal(0, sigma_T, n_points)

        try:
            poly = Chebyshev.fit(p_pert, T_pert, degree)
            deriv = poly.deriv()
            mu_samples[i, :] = deriv(p_meas)
        except Exception:
            mu_samples[i, :] = np.nan

    rel_err = (mu_samples - mu_true[np.newaxis, :]) / mu_true[np.newaxis, :] # relative error distribution at each point

    # drop extremities (first and last point)
    interior = slice(1, -1) if n_points > 2 else slice(None)
    rel_err_int = rel_err[:, interior]

    # clip outliers (|rel_err| > 5) that arise when mu_true is near zero
    # these extreme values destabilize the std estimate
    rel_err_int = np.clip(rel_err_int, -5, 5)

    # expanded uncertainty (k=1.96): 1.96 * std of the relative error
    point_rel_errs = k * np.nanstd(rel_err_int, axis=0)
    mean_rel_err = np.nanmean(point_rel_errs)

    return mean_rel_err, point_rel_errs


def _conv_uncert_for_slope(slope, dp, T_in, k=1.96):
    """
    Conventional relative expanded uncertainty for a given slope, pressure drop, and inlet temperature.
    T_out = T_in - slope * dp.
    """
    delta_T = slope * dp # absolute temperature change
    T_out = max(T_in - delta_T, 20.0) # outlet cannot be < 20 K
    return conventional_relative_uncertainty(
        mu_jt=slope, delta_T=delta_T, delta_p=dp,
        T_in=T_in, T_out=T_out)


def run_analysis():
    """
    Produce a figure of the theoretical JT measurement uncertainty
    as a function of the isenthalpic line slope (|dT/dp|).
    """
    plt.style.use(STYLE)
    p_in = 10.0   # MPa
    p_out = 1.5   # MPa
    dp = p_in - p_out  # 8.5 MPa

    slopes = np.concatenate([
        np.linspace(0.02, 0.5, 15),
        np.linspace(0.6, 2.0, 12),
        np.linspace(2.5, 8.0, 12),
        np.linspace(9.0, 15.0, 8),
    ])
    slopes_mc = np.concatenate([
        np.linspace(0.05, 0.5, 12),
        np.linspace(0.75, 2.0, 8),
        np.linspace(2.5, 8.0, 7),
        np.linspace(9.0, 15.0, 5),
    ])

    T_ref = 80.0  # representative inlet temperature

    # conventional uncertainty
    conv = {}
    for T in [65, 80, 160]:
        conv[T] = np.array([_conv_uncert_for_slope(s, dp, T) * 100
                            for s in slopes])

    # conventional curve on the MC slope grid
    conv_mc_grid = np.array([_conv_uncert_for_slope(s, dp, T_ref) * 100
                             for s in slopes_mc])

    # MC uncertainty
    n_points_list = [4, 5, 6, 8]
    curvatures = [0.0, 0.05]

    mc_results = {}
    for curv in curvatures:
        for n_pts in n_points_list:
            key = (n_pts, curv)
            errs = []
            for s in slopes_mc:
                me, _ = monte_carlo_uncertainty(
                    slope=s, p_in=p_in, p_out=p_out,
                    n_points=n_pts, T_in=T_ref,
                    curvature=curv, n_iter=2000)
                errs.append(me * 100)
            mc_results[key] = np.array(errs)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.5, 2.9), sharey=True)
    
    ax_a.fill_between(slopes, conv[65], conv[160],
                      color='0.75', alpha=0.35,
                      label=r'$T_{\rm in}$ = 65$\,$...$\,$160 K')
    ax_a.plot(slopes, conv[80], 'k-', lw=1.5, label=r'$T_{\rm in}$ = 80 K')

    # reference lines
    ax_a.axhline(5, color='0.4', ls=':', lw=0.8)
    ax_a.axhline(2, color='0.4', ls=':', lw=0.8)
    ax_a.text(13.5, 5.5, r'5 \%', fontsize=8.5, color='0.35', va='bottom')
    ax_a.text(13.5, 2.2, r'2 \%', fontsize=8.5, color='0.35', va='bottom')

    # fluid regions
    ax_a.axvspan(0, 0.5, alpha=0.08, color='C0')
    ax_a.axvspan(4, 8, alpha=0.06, color='C3')
    ax_a.text(0.55, 0.052, 'He', fontsize=8.5, color='C0', ha='center', weight='bold')
    ax_a.text(6.0, 0.052, r'N$_2$, Ar', fontsize=8.5, color='C3', ha='center', weight='bold')

    ax_a.set_xlabel(r'$|\mu_{\rm JT}|$/(K$\cdot$MPa$^{-1}$)')
    ax_a.set_ylabel(r'$U_{\rm r}(\mu_{\rm JT})$/\% ($k = 1.96$)')
    ax_a.set_title('(a) conventional')
    ax_a.set_yscale('log')
    ax_a.set_ylim(0.04, 100)
    ax_a.set_xlim(0, 15)
    ax_a.legend(fontsize=8.5, loc='upper right', frameon=True,
                edgecolor='black', facecolor='white', framealpha=1.0)

    colors = ['#4B0082', '#1E90FF', '#2E8B57', '#DAA520']
    ax_b.plot(slopes_mc, conv_mc_grid, color='0.55', lw=1.2, ls='-', label='conventional', zorder=1)

    # MC: solid = linear isenthalp, dashed = curved
    for ni, n_pts in enumerate(n_points_list):
        ax_b.plot(slopes_mc, mc_results[(n_pts, 0.0)],
                  '-', color=colors[ni], lw=1.4,
                  label=f'$n = {n_pts}$', zorder=2)
        ax_b.plot(slopes_mc, mc_results[(n_pts, 0.05)],
                  '--', color=colors[ni], lw=1.2, zorder=2)

    # reference lines
    ax_b.axhline(5, color='0.4', ls=':', lw=0.8)
    ax_b.axhline(2, color='0.4', ls=':', lw=0.8)
    ax_b.text(13.5, 5.5, r'5 \%', fontsize=8.5, color='0.35', va='bottom')
    ax_b.text(13.5, 2.2, r'2 \%', fontsize=8.5, color='0.35', va='bottom')

    # fluid regions
    ax_b.axvspan(0, 0.5, alpha=0.08, color='C0')
    ax_b.axvspan(4, 8, alpha=0.06, color='C3')
    ax_b.text(0.55, 0.052, 'He', fontsize=8.5, color='C0', ha='center', weight='bold')
    ax_b.text(6.0, 0.052, r'N$_2$, Ar', fontsize=8.5, color='C3',ha='center', weight='bold')

    # legends
    handles1, labels1 = ax_b.get_legend_handles_labels()
    style_handles = [Line2D([0], [0], color='0.4', ls='-', lw=1.3),
                     Line2D([0], [0], color='0.4', ls='--', lw=1.3)]
    leg1 = ax_b.legend(handles=style_handles, labels=['linear', 'curved'],
                fontsize=9, loc='upper right', frameon=True, edgecolor='black',
                facecolor='white', framealpha=1.0)
    ax_b.add_artist(leg1)
    ax_b.legend(handles1, labels1, fontsize=9, loc='upper center',
                       bbox_to_anchor=(0.45, 1.0), title='points / isenthalp',
                       title_fontsize=9, frameon=True, edgecolor='black',
                       facecolor='white', framealpha=1.0)

    ax_b.set_xlabel(r'$|\mu_{\rm JT}|$/(K$\cdot$MPa$^{-1}$)')
    ax_b.set_title('(b) Monte Carlo')
    ax_b.set_yscale('log')
    ax_b.set_xlim(0, 15)

    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(os.path.join(OUT_DIR, 'theoretical_uncertainty_combined.pdf'))


if __name__ == '__main__':
    run_analysis()
    plt.show()
