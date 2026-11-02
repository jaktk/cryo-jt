"""
Heat gained by the fluid between the inlet (TT101) and outlet (TT102)
thermometers, modelled for every point from the measured temperatures, and its
effect on mu_JT. The model does not use the EOS.

    Q = Q_rad + Q_out

  * Q_rad: radiation from the thermal shield onto the capillary,
        sigma * EPSILON * pi * D_CAP * L_CAP * mean(T_s^4 - T_c(z)^4),
    with T_s the highest shield temperature (TT009) and the capillary
    temperature T_c linear from TT101 to TT102 along its length. The capillary
    is treated as bare (EPSILON for unpolished stainless steel), so the
    multi-layer insulation is given no credit and Q_rad is an upper estimate.
  * Q_out: conduction down the outlet tube from the cryostat outlet (TT010) to
    the TT102 block, k * A_OUT * (TT010 - TT102) / L_OUT, with the outlet tube
    of the design calculation and the conductivity of 304 stainless steel
    (NIST cryogenic material properties).

A heat Q_i gained at point i warms the outlet by Q_i / (m_i c_p,i), with the
mass flow rate m_i from the capillary design calculation of the isenthalp at
the measured outlet pressure and c_p at the outlet state. Removing it and
refitting with the same degree and weights changes mu_JT by

    delta_Q = -d/dp [fit of Q_i / (m_i c_p,i)],

the local-derivative bias of mu_JT, which calculate_jt_coefficients.py treats
as the half-width of a rectangular distribution.
"""
import os
import numpy as np
import pandas as pd
from chebyshev_fit import fit_chebyshev
from repo_root import REPO_ROOT

MF_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'capillary_mass_flow')
SIGMA = 5.670374e-8                                   # W m^-2 K^-4
EPSILON = 0.1                                         # bare stainless steel
D_CAP, L_CAP = 1.0e-3, 6.0                            # m, capillary outer diameter and length
A_OUT = np.pi / 4 * ((4.0e-3)**2 - (3.0e-3)**2)       # m^2, outlet tube wall (4 mm x 3 mm)
L_OUT = 1.0                                           # m, outlet tube length
SS304_K = [-1.4087, 1.3982, 0.2543, -0.6260, 0.2334, 0.4256, -0.4658, 0.1650, -0.0199]

# Capillary design calculation of each measured isenthalp
GRIDS = {
    'Nitrogen_150K_6MPa': 'N2_D0.365_L6_p6.0_T149.8.csv',
    'Nitrogen_160K_5MPa': 'N2_D0.365_L6_p5.0_T160.0.csv',
    'Nitrogen_160K_6MPa': 'N2_D0.365_L6_p6.0_T160.0.csv',
    'Nitrogen_160K_7MPa': 'N2_D0.365_L6_p7.0_T160.0.csv',
    'Nitrogen_160K_9MPa': 'N2_D0.365_L6_p9.0_T160.0.csv',
    'Nitrogen_160K_12MPa': 'N2_D0.365_L6_p12.0_T160.0.csv',
    'Argon_180K_5MPa': 'Ar_D0.365_L6_p5.0_T180.0.csv',
    'Argon_180K_12MPa': 'Ar_D0.365_L6_p12.0_T180.0.csv',
    'Helium_65K_5MPa': 'He_D0.365_L6_p5.0_T65.0.csv',
    'Helium_65K_7MPa': 'He_D0.365_L6_p7.0_T65.0.csv',
    'Helium_141K_10MPa': 'He_D0.365_L6_p10.0_T140.6.csv',
    'Helium-Neon_65K_5MPa': 'HeNe_D0.365_L6_p5.0_T65.0_0.215He.csv',
    'Helium-Neon_65K_7MPa': 'HeNe_D0.365_L6_p7.0_T65.0_0.213He.csv',
    'Helium-Neon_65K_8MPa': 'HeNe_D0.365_L6_p8.0_T65.0_0.467He.csv',
    'Helium-Neon_65K_10MPa': 'HeNe_D0.365_L6_p10.0_T65.1_0.311He.csv',
    'Helium-Neon_66K_10MPa': 'HeNe_D0.365_L6_p10.0_T65.7_0.461He.csv',
    'Helium-Neon_80K_7MPa': 'HeNe_D0.365_L6_p7.0_T80.0_0.467He.csv',
    'Helium-Neon_80K_8MPa': 'HeNe_D0.365_L6_p8.0_T80.4_0.327He.csv',
    'Helium-Neon_80K_10MPa': 'HeNe_D0.365_L6_p10.0_T79.9_0.467He.csv',
    'Nitrogen-Helium_140K_8MPa': 'HeN2_D0.365_L6_p8.0_T140.0_0.503He.csv',
    'Nitrogen-Helium_160K_5MPa': 'HeN2_D0.365_L6_p5.0_T160.0_0.145He.csv',
    'Nitrogen-Helium_160K_8MPa': 'HeN2_D0.365_L6_p8.0_T160.0_0.145He.csv',
}


def ss304_conductivity(T):
    """ Thermal conductivity of 304 stainless steel in W/(m K), NIST fit (4 to 300 K). """
    x = np.log10(T)
    return 10**sum(c * x**i for i, c in enumerate(SS304_K))


def mass_flow(stem, p_out):
    """ Mass flow rate in kg/s at the outlet pressures p_out (MPa) from the design calculation. """
    grid = pd.read_csv(os.path.join(MF_DIR, GRIDS[stem])).sort_values('p2/MPa')
    return np.interp(p_out, grid['p2/MPa'], grid['m/g_per_s']) / 1e3


def heat_gain(df):
    """ Radiative and outlet-conduction heat gained by the fluid at each point, in W. """
    t_in, t_out = df['TT101/K'].to_numpy(float), df['TT102/K'].to_numpy(float)
    t_shield, t_exit = df['TT009/K'].to_numpy(float), df['TT010/K'].to_numpy(float)
    z = np.linspace(0.0, 1.0, 51)
    t_cap = t_in[:, None] + (t_out - t_in)[:, None] * z[None, :]
    q_rad = SIGMA * EPSILON * np.pi * D_CAP * L_CAP * np.mean(t_shield[:, None]**4 - t_cap**4, axis=1)
    k = 0.5 * (ss304_conductivity(t_exit) + ss304_conductivity(t_out))
    q_out = k * A_OUT * (t_exit - t_out) / L_OUT
    return q_rad, q_out


def heat_bias(stem, df, fluid_props, degree, weights):
    """ Local-derivative bias of mu_JT (K/MPa) at every point from the modelled
    heat gain: the change of the fitted slope when the outlet warming
    Q_i / (m_i c_p,i) is removed. Returns the bias and the heat gain in W. """
    p, t = df['PT102/MPa'].to_numpy(float), df['TT102/K'].to_numpy(float)
    q_rad, q_out = heat_gain(df)
    q = q_rad + q_out
    cp = np.array([fluid_props.get_cp(a, b) for a, b in zip(p, t)])
    warming = q / (mass_flow(stem, p) * cp)
    bias = -fit_chebyshev(p, warming, degree, weights).deriv()(p)
    return bias, q
