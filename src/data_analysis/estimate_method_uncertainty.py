"""
Estimate the relative method uncertainty of mu_JT from the pure-fluid validation.

The N2 and Ar points scatter around their reference equations of state
(Span et al. 2000; Tegeler et al. 1999) by more than the sensor-level Monte
Carlo uncertainty. The excess is attributed to systematic effects that the
Monte Carlo does not model (flow-dependent heat exchange, inlet-thermometer
offset, differentiation of few-point isenthalps). It is estimated as

    u_method^2 = RMS(dev)^2 - mean(u_MC^2 + u_sys^2 + u_EOS^2)

with all terms relative and at k = 1. Neither EOS publication states an
uncertainty for mu_JT; u_EOS is taken from their isobaric heat capacity
tolerance diagrams (Span et al. Fig. 41, Tegeler et al. Fig. 44), since
mu_JT = (T (dv/dT)_p - v) / c_p.

Run after calculate_jt_coefficients.py; the result is the constant
U_METHOD_REL used there.
"""
import glob
import os
import numpy as np
import pandas as pd
from get_git_root import get_git_root
from calculate_jt_coefficients import K_TARGET, systematic_rel_std_uncertainty


def eos_rel_std_uncertainty(fluid, p_mpa):
    """ Relative standard uncertainty (k = 1) of c_p from the EOS tolerance diagrams (expanded values / 2). """
    if fluid == 'Nitrogen':
        return 0.8 / 100 / 2   # region C, (100 to 160) K below 30 MPa
    if fluid == 'Argon':
        return (2.0 if p_mpa >= 5.5 else 0.3) / 100 / 2   # region C above about 5.5 MPa at (155 to 175) K, else A
    raise ValueError(fluid)


def main():
    root = get_git_root(os.getcwd())
    rows = []
    for fluid in ('Nitrogen', 'Argon'):
        for f in sorted(glob.glob(os.path.join(root, 'data', 'derived_data', 'jt_coeffs', f'{fluid}_*_JT.csv'))):
            df = pd.read_csv(f)
            df = df[df['used'] == True]
            mu_eos = df['JT_eos/(K/MPa)'].values
            mu = df['JT_meas/(K/MPa)'].values
            dev = (mu - mu_eos) / mu_eos
            u_mc = df['JT_UNC_MC/(K/MPa)'].values / K_TARGET / np.abs(mu)
            u_sys = systematic_rel_std_uncertainty(fluid)
            u_eos = np.array([eos_rel_std_uncertainty(fluid, p) for p in df['p/MPa'].values])
            rows.append(pd.DataFrame({'dev': dev, 'u_mc': u_mc, 'u_sys': u_sys, 'u_eos': u_eos}))
    d = pd.concat(rows)
    rms = np.sqrt(np.mean(d['dev']**2))
    explained = np.mean(d['u_mc']**2 + d['u_sys']**2 + d['u_eos']**2)
    u_method = np.sqrt(rms**2 - explained)
    print(f"N2 + Ar points: {len(d)}")
    print(f"mean relative deviation: {100 * d['dev'].mean():+.2f} %")
    print(f"RMS relative deviation:  {100 * rms:.2f} %")
    print(f"RMS of modelled terms:   {100 * np.sqrt(explained):.2f} % "
          f"(MC {100 * np.sqrt(np.mean(d['u_mc']**2)):.2f}, bias {100 * np.sqrt(np.mean(d['u_sys']**2)):.2f}, "
          f"EOS {100 * np.sqrt(np.mean(d['u_eos']**2)):.2f})")
    print(f"u_method (k = 1):        {100 * u_method:.2f} %")


if __name__ == '__main__':
    main()
