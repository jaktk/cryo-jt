"""
Measurement uncertainties of the outlet pressure and temperature, shared by
every script that needs them. No REFPROP needed.

All values are expanded with K_TARGET (approx. 95 %) unless stated otherwise.

Temperature (outlet thermometer TT102, Cernox X93303), three independent
contributions combined in quadrature:
  * the Cernox calibration: the CX-1050 column of the uncertainty table of the
    Lake Shore calibration report of X93303 (data/cernox_calibration_data/
    X93303/X93303.pdf, k = 2), interpolated linearly in temperature;
  * the acquisition module CABTR: a fit to the specifications of its manual;
  * the calibration polynomial T(R): 2 sqrt(N/(N-n)) DT_rms per band of the
    Lake Shore fit.
Pressure: Mensor CPT 6100, accuracy 0.01 % of the 13.7 MPa full span, stated
by the manufacturer at k = 2.

The per-point uncertainty combines in quadrature the sensor-chain value and
the statistical uncertainty of the mean over the steady-state window
(the *_EXP_UNC columns, 2 s / sqrt(N) with N = N_WINDOW samples). The
statistical standard uncertainty is scaled by t / 2, with t the 0.975 quantile
of Student's t for N - 1 degrees of freedom, so that its expanded value
corresponds to a level of confidence of about 95 % for the short window.
"""
import numpy as np
from scipy import stats

K_TARGET = 2.0   # coverage factor used throughout the analysis
K_SOURCE = 2.0   # coverage factor of the sensor and calibration specifications
N_WINDOW = 5     # samples averaged per steady-state point (get_pT_pairs.py)
STUDENT_SCALE = stats.t.ppf(0.975, N_WINDOW - 1) / 2.0

P_FULL_SCALE = 13.7  # MPa
U_P = 0.0001 * P_FULL_SCALE / K_SOURCE  # MPa, standard (k = 1) pressure uncertainty

# Calibration report of X93303, uncertainty table, CX-1050 column (k = 2)
CERNOX_T = np.array([1.4, 4.2, 10.0, 20.0, 30.0, 50.0, 100.0, 300.0, 400.0])   # K
CERNOX_U = np.array([4.0, 4.0, 4.0, 8.0, 9.0, 12.0, 17.0, 46.0, 74.0])         # mK

# Lake Shore Chebyshev fit of X93303 per temperature band
POLY_BANDS = [
    {"Tmin": 1.4, "Tmax": 14.1, "N": 31, "n": 9, "DTrms": 0.93},
    {"Tmin": 14.1, "Tmax": 80.0, "N": 31, "n": 7, "DTrms": 1.60},
    {"Tmin": 80.0, "Tmax": 325.0, "N": 32, "n": 8, "DTrms": 5.84},
]


def cernox_uncertainty(T):
    """ Cernox calibration uncertainty in mK (k = K_SOURCE). """
    return np.interp(T, CERNOX_T, CERNOX_U)


def cabtr_uncertainty(T):
    """ CABTR acquisition-module uncertainty in mK (k = K_SOURCE), fit to the
    specifications of the user manual, 20 K to 300 K. """
    return 2e-05 * np.asarray(T)**2 + 0.0695 * np.asarray(T) - 0.1001


def polynomial_fit_uncertainty(T):
    """ Uncertainty of the calibration polynomial T(R) in mK (k = K_SOURCE). """
    def band(t):
        for b in POLY_BANDS:
            if b["Tmin"] <= t < b["Tmax"]:
                return K_SOURCE * np.sqrt(b["N"] / (b["N"] - b["n"])) * b["DTrms"]
        return 0.0
    return np.vectorize(band)(T)


def temperature_chain_uncertainty(T):
    """ Sensor-chain uncertainty of TT102 in K, expanded with K_TARGET. """
    T = np.asarray(T, dtype=float)
    u = np.sqrt(cernox_uncertainty(T)**2 + cabtr_uncertainty(T)**2 + polynomial_fit_uncertainty(T)**2)
    return u * 1e-3 * (K_TARGET / K_SOURCE)


def pressure_chain_uncertainty():
    """ Sensor uncertainty of PT102 in MPa, expanded with K_TARGET. """
    return K_TARGET * U_P


def combine(statistical_exp_unc, chain_exp_unc):
    """ Per-point expanded uncertainty (k = K_TARGET) from the window
    statistics (2 s / sqrt(N)) and the sensor chain, in quadrature. """
    u_stat = STUDENT_SCALE * np.asarray(statistical_exp_unc, dtype=float) / 2.0
    u_chain = np.asarray(chain_exp_unc, dtype=float) / K_TARGET
    return K_TARGET * np.sqrt(u_stat**2 + u_chain**2)


def point_uncertainties(df):
    """ Expanded (k = K_TARGET) uncertainties of the outlet pressure (MPa) and
    temperature (K) of every point of a p_T_pairs table. """
    T = df['TT102/K'].to_numpy(dtype=float)
    p_unc = combine(df['PT102/MPa_EXP_UNC'].to_numpy(dtype=float), pressure_chain_uncertainty())
    t_unc = combine(df['TT102/K_EXP_UNC'].to_numpy(dtype=float), temperature_chain_uncertainty(T))
    return p_unc, t_unc
