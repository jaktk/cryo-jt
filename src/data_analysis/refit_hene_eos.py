"""
Refit of the He-Ne mixture model to the measured Joule-Thomson coefficients.

The binary parameters of the He-Ne model in REFPROP/FLUIDS/HMX.BNC (model BHN:
the reducing parameters betaT, gammaT, betaV, gammaV and the weight F of the
departure function) are refitted to the measured mu_JT of the He-Ne points
used in the paper, weighted by their standard uncertainty U/2. The shape of the
departure function is kept. Two variants are fitted, each to all He-Ne points
and to the well-conditioned isenthalps only (N_s >= 4 and chi-square test
passed):

  * free: the five parameters within BOUNDS;
  * phase envelope held: bubble points of the published model on a grid of
    temperatures and liquid compositions are added as pseudo-data, with
    tolerances TOL_P (relative pressure) and TOL_Y (vapour mole fraction), as
    a stand-in for the phase-equilibrium data to which the model was fitted.

For each fit the script reports the AAD (referred to the measured value, as
in the paper) and the number of points outside U, and the largest change from
the published model in the bubble-point pressure, the vapour composition, and
the density at the measured states. The parameters are changed in memory only
(REFPROP SETKTVdll); no file of the REFPROP installation is modified.

Output: data/derived_data/jt_coeffs/hene_eos_refit.csv

Run after calculate_jt_coefficients.py.
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from FluidProps import FluidProps
from repo_root import REPO_ROOT

JT_DIR = os.path.join(REPO_ROOT, 'data', 'derived_data', 'jt_coeffs')
VLE_T = [30.0, 33.0, 36.0, 39.0, 42.0]   # K, bubble-point grid of the constraint
VLE_X = [0.005, 0.01, 0.02]              # liquid helium mole fraction
TOL_P, TOL_Y = 0.01, 0.005
BOUNDS = ([0.5, 0.5, 0.5, 0.5, 0.0], [1.5, 1.5, 1.5, 1.5, 3.0])


def load_points():
    """ He-Ne points used for mu_JT, with the isenthalp-mean composition used
    for the EOS reference in calculate_jt_coefficients.py. """
    summary = pd.read_csv(os.path.join(JT_DIR, 'jt_coefficients_summary.csv'))
    summary = summary[summary['fluid'] == 'Helium-Neon']
    rows = []
    for _, meta in summary.iterrows():
        df = pd.read_csv(os.path.join(JT_DIR, meta['filename'].replace('.csv', '_JT.csv')))
        df = df[df['used'] == True]
        well_conditioned = meta['n_points'] >= 4 and meta['chi2_pass'] == True
        for _, r in df.iterrows():
            rows.append({'x1': meta['x1'], 'p': r['p/MPa'], 'T': r['T/K'],
                         'mu': r['JT_meas/(K/MPa)'], 'U': r['JT_UNC/(K/MPa)'],
                         'mu_eos': r['JT_eos/(K/MPa)'], 'well_conditioned': well_conditioned})
    return pd.DataFrame(rows)


class HeNeModel:
    def __init__(self, points):
        self.points = points
        self.fp = FluidProps(('Helium', 'Neon'))
        self.fp.set_composition_from_1st_fraction(0.3)
        self.fp.get_JT_coefficient(5.0, 60.0)   # loads the mixture and its BHN parameters
        ktv = self.fp.RP.GETKTVdll(1, 2)
        self.hmodij, self.hfmix = ktv.hmodij, ktv.hFmix
        self.published = np.array(ktv.fij[:5])

    def set_parameters(self, f):
        out = self.fp.RP.SETKTVdll(1, 2, self.hmodij, list(f) + [0.0], self.hfmix)
        if out.ierr != 0:
            raise RuntimeError(out.herr)

    def _at_points(self, prop):
        out = []
        for r in self.points.itertuples():
            self.fp.set_composition_from_1st_fraction(r.x1)
            out.append(self.fp.RP.REFPROPdll(self.fp.fluid, 'PT', prop, self.fp.MASS_BASE_SI, 0, 0,
                                             r.p * 1e6, r.T, self.fp.z).Output[0])
        return np.array(out)

    def mu(self):
        return self._at_points('JT') * 1e6   # K/MPa

    def density(self):
        return self._at_points('D')

    def bubble_points(self):
        p, y = [], []
        for T in VLE_T:
            for x in VLE_X:
                s = self.fp.RP.SATTdll(T, [x, 1 - x], 1)
                p.append(s.P if s.ierr == 0 else np.nan)
                y.append(s.y[0] if s.ierr == 0 else np.nan)
        return np.array(p), np.array(y)


def main():
    points = load_points()
    model = HeNeModel(points)
    model.set_parameters(model.published)
    mu0, rho0 = model.mu(), model.density()
    pb0, yb0 = model.bubble_points()
    if not np.allclose(mu0, points['mu_eos'], rtol=1e-6):
        raise RuntimeError('EOS mu_JT differs from the pipeline output')
    ok = np.isfinite(pb0)

    def residuals(f, subset, hold_envelope):
        model.set_parameters(f)
        r = ((model.mu() - points['mu']) / (points['U'] / 2)).values[subset]
        if hold_envelope:
            pb, yb = model.bubble_points()
            rp = np.where(np.isfinite(pb), (pb - pb0) / pb0 / TOL_P, 10.0)[ok]
            ry = np.where(np.isfinite(yb), (yb - yb0) / TOL_Y, 10.0)[ok]
            r = np.concatenate([r, rp, ry])
        return r

    def row(label, f):
        model.set_parameters(f)
        mu = model.mu()
        pb, yb = model.bubble_points()
        dev = np.abs(mu - points['mu']) / np.abs(points['mu'])
        outside = np.abs(mu - points['mu']) > points['U']
        wc = points['well_conditioned'].values
        return {'fit': label,
                'betaT': f[0], 'gammaT': f[1], 'betaV': f[2], 'gammaV': f[3], 'F': f[4],
                'AAD_all_perc': 100 * dev.mean(), 'outside_U_all': int(outside.sum()), 'n_all': len(points),
                'AAD_wc_perc': 100 * dev[wc].mean(), 'outside_U_wc': int(outside[wc].sum()), 'n_wc': int(wc.sum()),
                'max_dp_bubble_perc': 100 * np.nanmax(np.abs(pb - pb0) / pb0),
                'max_dy_bubble': np.nanmax(np.abs(yb - yb0)),
                'max_drho_perc': 100 * np.max(np.abs(model.density() - rho0) / rho0)}

    rows = [row('published', model.published)]
    everything = np.ones(len(points), dtype=bool)
    wc = points['well_conditioned'].values
    for label, subset, hold in [('all, free', everything, False),
                                ('all, phase envelope held', everything, True),
                                ('well-conditioned, free', wc, False),
                                ('well-conditioned, phase envelope held', wc, True)]:
        res = least_squares(residuals, model.published, args=(subset, hold),
                            x_scale=0.05, bounds=BOUNDS)
        rows.append(row(label, res.x))
    model.set_parameters(model.published)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(JT_DIR, 'hene_eos_refit.csv'), index=False)
    print(out.round(4).to_string(index=False))


if __name__ == '__main__':
    main()
