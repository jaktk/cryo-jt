import os
import json
import numpy as np
import pandas as pd
from scipy import optimize, interpolate
from numpy.polynomial.chebyshev import Chebyshev
from sklearn.metrics import r2_score
from FluidProps import FluidProps
from get_git_root import get_git_root

# Coverage factors. The sensor, acquisition, and calibration specifications are
# expanded at K_SOURCE; the analysis reports everything at K_TARGET (95 %).
K_TARGET = 1.96   # coverage factor used throughout the analysis (95 % level of confidence)
K_SOURCE = 2.0    # coverage factor of the sensor / calibration source specifications
N_MC_SAMPLES = 2000  # Monte Carlo realisations for the per-point uncertainty


def chebyshev_fit_weights(pressures, temperatures, p_uncertainties, t_uncertainties):
    """
    Weights for the Chebyshev least-squares fit. The pressure uncertainty is
    folded into an equivalent temperature uncertainty through the local slope
    dT/dp, so that weight = 1 / (u_T^2 + (u_p * dT/dp)^2). Returns normalised
    weights (sum to one)."""
    order = np.argsort(pressures)
    slope = np.empty_like(temperatures, dtype=float)
    slope[order] = np.gradient(temperatures[order], pressures[order])
    weights = 1.0 / (t_uncertainties**2 + (p_uncertainties * slope)**2)
    return weights / np.sum(weights)


class CernoxCal(object):
    serNum = ''

    def __init__(self, serNum):
        self.serNum = serNum
        self.coefficients = self._set_coefficients()
        self.tbl = self._set_tbl()
        self.dat = self._set_dat()

    def _set_coefficients(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.cof')
        c_dict = {}
        with open(fpath, 'r') as fh:
            coefs = []
            limits = ''
            while True:
                line = fh.readline()
                if not line:
                    c_dict[limits] = {'Zlower': Zlower, 'Zupper': Zupper, 'coefs': coefs}
                    break
                l = line.strip().split()
                if l[0] == 'Lower' and l[1] == 'Resist.':
                    limits = l[-1]
                elif l[0] == 'Upper' and l[1] == 'Resist.':
                    limits = limits + ',' + l[-1]
                elif l[0] == 'Zlower':
                    Zlower = float(l[-1])
                elif l[0] == 'Zupper':
                    Zupper = float(l[-1])
                elif list(l[0])[0] == 'C' and list(l[0])[1] == '(':
                    coefs.append(float(l[-1]))
                elif coefs and list(l[0])[0] != 'C':
                    c_dict[limits] = {'Zlower': Zlower, 'Zupper': Zupper, 'coefs': coefs}
                    coefs = []
        return c_dict

    def _set_tbl(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.tbl')
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def _set_dat(self):
        fpath = os.path.join(get_git_root(os.getcwd()),
                             'data',
                             'cernox_calibration_data',
                             self.serNum,
                             f'{self.serNum}.dat')
        return pd.read_csv(fpath, sep=r'\s+', header=[0,1])

    def get_T_from_coefs(self, R):
        T = 0
        for key, value in self.coefficients.items():
            L, U = key.split(',')
            L, U = float(L), float(U)
            if R >= L and R < U:
                k = ((np.log10(R) - value['Zlower']) - (value['Zupper'] - np.log10(R))) / (value['Zupper'] - value['Zlower'])
                for i, c in enumerate(value['coefs']):
                    T += c * np.cos(i * np.arccos(k))
        return T

    def get_T_from_tbl(self, R):
        f = interpolate.interp1d(self.tbl[('Resistance', '(Ohms)')], self.tbl[('Temp.', '(Kelvin)')])
        return f(R)

    def get_interp_abs_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Ti = self.get_T_from_tbl(R)
        return Tc - Ti

    def get_interp_rel_error(self, R):
        Tc = self.get_T_from_coefs(R)
        Ti = self.get_T_from_tbl(R)
        return (Tc - Ti) / Ti

    def get_deltaT_dat_cof(self):
        T, DT = [], []
        for row in self.dat.iterrows():
            Tdat = row[1][('Temperature','(Kelvin)')]
            T.append(Tdat)
            DT.append(Tdat - self.get_T_from_coefs(row[1][('Resistance','(Ohms)')]))
        return T, DT


class TempUncertainty(object):
    sensor = ''
    def __init__(self, sensor):
        assert sensor in ['X93303','X115143']
        self.sensor = sensor
        self.cernox = CernoxCal(sensor)

    def _get_R_from_T(self, T):
        dfx0 = self.cernox.dat.iloc[(self.cernox.dat[('Temperature','(Kelvin)')]-T).abs().argsort()[:1]]
        res = optimize.minimize(fun = lambda x: np.abs(self.cernox.get_T_from_coefs(x) - T),
                                x0 = dfx0[('Resistance','(Ohms)')])
        return res.x[0]

    def get_cernox_uncertainty(self, T):
        """ Uncertainty in mK for temperature from 20 K to 300 K """
        return -5e-05 * T**2 + 0.1281 * T + 6.4603

    def get_cabtr_uncertainty(self, T):
        """
        Uncertainty in mK for temperature from 20 K to 300 K.
        Fit to data from CABTR user manual.
        """
        return 2e-05 * T**2 + 0.0695 * T - 0.1001

    def get_polynomial_fit_uncert(self, T):
        """ uncertainty in mK from fitting Chebyshev polynomial """
        d = {
        'X93303': [{'Tmin': 1.4 , 'Tmax': 14.1 , 'N': 31, 'n': 9, 'DTrms': 0.93},
                  {'Tmin': 14.1, 'Tmax': 80.0 , 'N': 31, 'n': 7, 'DTrms': 1.60},
                  {'Tmin': 80.0, 'Tmax': 325.0, 'N': 32, 'n': 8, 'DTrms': 5.84}],
        'X115143': [{'Tmin': 20.0, 'Tmax': 95.0, 'N': 28, 'n': 7, 'DTrms': 0.91},
                   {'Tmin': 95.0, 'Tmax': 325.0, 'N': 29, 'n': 9, 'DTrms': 3.45}],
        }

        for i in d[self.sensor]:
            if T >= i['Tmin'] and T < i['Tmax']:
                sigm2 = i['N'] / (i['N'] - i['n']) * i['DTrms']**2
                break
            else:
                sigm2 = 0
        return K_SOURCE * sigm2**0.5

    def __call__(self, T):
        """ return combined temperature uncertainty in K (expanded, k = K_TARGET, 95 %) """
        cabtr = self.get_cabtr_uncertainty(T)
        cernox = self.get_cernox_uncertainty(T)
        poly = self.get_polynomial_fit_uncert(T)
        # The sensor, acquisition, and calibration-polynomial contributions above
        # are expanded at k = K_SOURCE; sum them linearly and rescale the combined
        # chain to k = K_TARGET (95 % level of confidence) used across the analysis.
        return (cabtr + cernox + poly) * 1e-3 * (K_TARGET / K_SOURCE)


class JTCoefficientCalculator(object):
    def __init__(self):
        self.git_root = get_git_root(os.getcwd())
        self.data_dir = os.path.join(self.git_root,
                                     'data',
                                     'derived_data',
                                     'p_T_pairs')
        self.metadata_file = os.path.join(self.git_root,
                                          'data',
                                          'derived_data',
                                          'p_T_pairs',
                                          'index.json')
        self.output_dir = os.path.join(self.git_root,
                                       'data',
                                       'derived_data',
                                       'jt_coeffs')
        os.makedirs(self.output_dir, exist_ok=True)
        self.metadata = self._load_metadata()
        
    def _load_metadata(self):
        with open(self.metadata_file, 'r') as f:
            metadata = json.load(f)
        metadata_dict = {}
        for item in metadata:
            metadata_dict[item['filename']] = item
        return metadata_dict
    
    def _parse_fluid_info(self, filename):
        metadata = self.metadata.get(filename, {})
        fluid_name = metadata.get('fluid', '')
        fluid_type = fluid_name.split("-")
        if len(fluid_type) > 1:
            fluid_type = tuple(fluid_type)
            mixture = True
        else:
            fluid_type = fluid_type[0]
            mixture = False
        return fluid_type, mixture
    
    def _fit_polynomial_weighted(self,
                                 pressures,
                                 temperatures,
                                 p_uncertainties,
                                 t_uncertainties,
                                 max_degree=6):
        """
        Fit optimal Chebyshev polynomial to p-T data with weighted least squares.
        Uses measurement uncertainties to weight the fitting process.
        """
        n_points = len(pressures)
        max_degree = min(max_degree, n_points - 1)

        # Weights from the measurement uncertainties (pressure folded into an
        # equivalent temperature uncertainty through the local slope dT/dp).
        weights = chebyshev_fit_weights(pressures, temperatures,
                                        p_uncertainties, t_uncertainties)

        best_degree = 1
        best_r2 = -np.inf
        best_fit = None
        
        for degree in range(1, max_degree + 1):
            try:
                # Fit weighted Chebyshev polynomial
                poly_fit = Chebyshev.fit(pressures, temperatures, degree, w=weights)
                poly = Chebyshev(poly_fit.convert().coef)
                
                # Calculate R^2
                fitted_temps = poly(pressures)
                r2 = r2_score(temperatures, fitted_temps)
                
                # Store best fit
                if r2 > best_r2:
                    best_r2 = r2
                    best_degree = degree
                    best_fit = poly
        
                # Break if R^2 > 0.999
                if r2 > 0.999:
                    break
                    
            except Exception as e:
                print(f"Error fitting degree {degree}: {e}")
                continue
        
        return best_fit, best_degree, best_r2
    
    def _calculate_jt_coefficients(self, pressures, polynomial):
        """ Calculate JT coefficients at measurement points with polynomial derivative """
        poly_deriv = polynomial.deriv(m=1) # Get derivative of polynomial (this gives dT/dp)
        return poly_deriv(pressures)
    
    def _calculate_theoretical_jt(self,
                                  pressures,
                                  temperatures,
                                  fluid_props,
                                  composition=None):
        """ Calculate theoretical JT coefficients with EOS and REFPROP """
        jt_theoretical = np.zeros_like(pressures)
        
        for i, (p, t) in enumerate(zip(pressures, temperatures)):
            try:
                if composition is not None:
                    fluid_props.set_composition_from_1st_fraction(composition) # Set composition for mixtures
                jt_theoretical[i] = fluid_props.get_JT_coefficient(p, t)
                
            except Exception as e:
                print(f"Error calculating theoretical JT for p={p}, T={t}: {e}")
                jt_theoretical[i] = np.nan
        
        return jt_theoretical
    
    def _monte_carlo_uncertainty(self,
                                 pressures,
                                 temperatures,
                                 p_uncertainties,
                                 t_uncertainties,
                                 polynomial,
                                 n_samples=N_MC_SAMPLES,
                                 seed=20210101):
        """
        Monte Carlo expanded uncertainty (k = K_TARGET, 95 %) of the derived
        Joule-Thomson coefficient.

        Each measured point is perturbed by a normal (Gaussian) error whose
        standard deviation is the *standard* (k = 1) measurement uncertainty,
        i.e. the expanded uncertainty divided by K_TARGET. The perturbed data
        are refitted with the same weighted Chebyshev fit used for the central
        value, and the coefficient is re-evaluated by differentiation. The
        expanded uncertainty is K_TARGET times the standard deviation of the
        resulting distribution. A Gaussian (rather than uniform) error reflects
        the calibration-derived nature of the measurement uncertainties.
        """
        try:
            # Convert the expanded (k = K_TARGET) uncertainties to standard
            # (k = 1) deviations for sampling.
            p_sigma = p_uncertainties / K_TARGET
            t_sigma = t_uncertainties / K_TARGET
            degree = polynomial.degree()
            # Weight the refit as the central fit does, but only when the fit is
            # over-determined: for an exact interpolation (degree == n - 1) the
            # weights cannot change the interpolant and only destabilise the
            # weighted Vandermonde, so an unweighted refit is used there.
            weights = (chebyshev_fit_weights(pressures, temperatures,
                                             p_uncertainties, t_uncertainties)
                       if degree < len(pressures) - 1 else None)

            jt_samples = []
            rng = np.random.default_rng(seed)

            for _ in range(n_samples):
                # Gaussian perturbation within the standard measurement uncertainty
                p_sample = rng.normal(pressures, p_sigma)
                t_sample = rng.normal(temperatures, t_sigma)

                # Refit with the same weighted Chebyshev fit as the central value
                try:
                    poly_fit = Chebyshev.fit(p_sample, t_sample,
                                             degree, w=weights)
                    poly_sample = Chebyshev(poly_fit.convert().coef)

                    # Calculate JT coefficients
                    poly_deriv = poly_sample.deriv(m=1)
                    jt_sample = poly_deriv(pressures)

                    jt_samples.append(jt_sample)

                except Exception:
                    continue

            if jt_samples:
                jt_samples = np.array(jt_samples)
                # Expanded uncertainty: k = K_TARGET times the sample std (k = 1).
                jt_mc_uncertainty = K_TARGET * np.std(jt_samples, axis=0, ddof=1)
                return jt_mc_uncertainty
            else:
                return np.zeros_like(pressures)

        except Exception as e:
            print(f"Error in Monte Carlo uncertainty: {e}")
            return np.zeros_like(pressures)
    
    def process_file(self, filename):
        """ Process a single p_T_pairs file to calculate JT coefficients """
        filepath = os.path.join(self.data_dir, filename)
        
        print(f"Processing {filename}...")
        
        try:
            df = pd.read_csv(filepath) # Read data
            metadata = self.metadata.get(filename, {}) # Get metadata
            fluid_type, is_mixture = self._parse_fluid_info(filename) # Parse fluid information from metadata
            pressures = df['PT102/MPa'].values # Downstream pressure
            temperatures = df['TT102/K'].values # Downstream temperature
            
            # Per-point measurement uncertainty: max of the steady-state statistical
            # uncertainty (from the averaging window) and the sensor-chain expanded
            # uncertainty (k=1.96) at the point's pressure / temperature.
            p_uncertainties = np.array([max(p, 1.96 * 0.01/100 * 13.7) for p in df['PT102/MPa_EXP_UNC'].values])
            comb_T_unc = TempUncertainty("X93303")
            t_uncertainties = np.array([
                max(t_stat, comb_T_unc(T))
                for t_stat, T in zip(df['TT102/K_EXP_UNC'].values, temperatures)
            ])
            
            # Skip if insufficient data
            if len(pressures) < 3:
                print(f"  Insufficient data points: {len(pressures)}")
                return None
            
            fluid_props = FluidProps(fluid_type)
            
            # Get composition, if mixture
            composition = None
            x1_per_point = None
            if is_mixture:
                if 'x1' in df.columns and not df['x1'].isna().all():
                    x1 = df['x1'].values.astype(float)
                    x2 = (df['x2'].values.astype(float)
                          if 'x2' in df.columns else 1.0 - x1)
                    # Reject out-of-range gas-analyzer (GA244) readings
                    valid = np.abs(x1 + x2 - 1.0) < 0.02
                    if not valid.any():
                        valid = np.ones_like(x1, dtype=bool)
                    composition = x1[valid].mean()
                    n_rejected = int((~valid).sum())
                    if n_rejected:
                        print(f"  Rejected {n_rejected} out-of-range GA reading(s); "
                              f"x1 = {composition:.4f} (was {x1.mean():.4f} with outliers)")
                    # keep valid readings, assign the isenthalp mean to rejected points
                    x1_per_point = np.where(valid, x1, composition)
                    fluid_props.set_composition_from_1st_fraction(composition)
                else:
                    print(f"  Warning: No composition data found for mixture {filename}")
            
            # Fit polynomial to data with uncertainties
            poly, degree, r2 = self._fit_polynomial_weighted(pressures, temperatures, p_uncertainties, t_uncertainties)
            
            if poly is None:
                print(f"  Failed to fit polynomial")
                return None
            
            # Calculate JT coefficients with full uncertainty propagation
            jt_measured = self._calculate_jt_coefficients(pressures, poly)
            
            # Calculate Monte Carlo uncertainty
            jt_uncertainty = self._monte_carlo_uncertainty(
                pressures, temperatures, p_uncertainties, t_uncertainties, poly
            )
            
            # Calculate theoretical JT coefficients
            jt_theoretical = self._calculate_theoretical_jt(pressures, temperatures, fluid_props, composition)
            
            # Calculate relative errors of EOS, assuming the ground truth is: EOS for pure fluids & measurements for mixtures
            denominator = jt_measured if is_mixture else jt_theoretical
            relative_abs_error = np.where(jt_theoretical != 0,
                                          np.abs((jt_measured - jt_theoretical) / denominator * 100),
                                          np.nan)

            # Exclude the isenthalp extremities
            interior_mask = np.ones(len(pressures), dtype=bool)
            if len(pressures) > 2:
                order = np.argsort(pressures)
                interior_mask[order[0]] = False
                interior_mask[order[-1]] = False
            rel_err_interior = relative_abs_error[interior_mask]
            unc_interior = (jt_uncertainty / np.abs(jt_measured) * 100)[interior_mask]

            # Create results dictionary
            results = {
                'filename': filename,
                'fluid': fluid_type if isinstance(fluid_type, str) else f"{fluid_type[0]}-{fluid_type[1]}",
                'is_mixture': is_mixture,
                'x1': composition,
                'x2': None if not is_mixture else (1 - composition),
                'n_points': len(pressures),
                'polynomial_degree': degree,
                'r2_score': r2,
                'TT101_mean/K': metadata.get('mean_TT101_K', np.nan),
                'p/MPa': pressures,
                'T/K': temperatures,
                'x1_per_point': x1_per_point,
                'p_UNC/MPa': p_uncertainties,
                'T_UNC/K': t_uncertainties,
                'JT_meas/(K/MPa)': jt_measured,
                'JT_eos/(K/MPa)': jt_theoretical,
                'JT_UNC/(K/MPa)': jt_uncertainty,
                'used': interior_mask,
                'rel_abs_err_perc': relative_abs_error,
                'mean_rel_abs_err_perc': np.nanmean(rel_err_interior),
                'std_rel_abs_err_perc': np.nanstd(rel_err_interior),
                'rms_rel_abs_err_perc': np.sqrt(np.nanmean(rel_err_interior**2)),
                'mean_unc_perc': np.nanmean(unc_interior)
            }
            
            print(f"  Success: {len(pressures)} points, degree={degree}, R²={r2:.6f}")
            print(f"  Mean relative error (interior): {np.nanmean(rel_err_interior):.2f}%")
            print(f"  Mean uncertainty (interior): {np.nanmean(unc_interior):.2f}%")
            
            return results
            
        except Exception as e:
            print(f"  Error processing {filename}: {e}")
            return None
    
    def process_all_files(self):
        """ Process all files in the p_T_pairs directory """
        # Get all CSV files from metadata
        csv_files = list(self.metadata.keys())
        csv_files.sort()
        
        print(f"Found {len(csv_files)} files in metadata")
        
        all_results = []
        
        for csv_file in csv_files:
            # Check if file exists
            filepath = os.path.join(self.data_dir, csv_file)
            if not os.path.exists(filepath):
                print(f"Warning: File {csv_file} listed in metadata but not found in directory")
                continue
                
            result = self.process_file(csv_file)
            if result is not None:
                all_results.append(result)
        
        return all_results
    
    def save_results(self, results):
        """ Save results to files """
        if not results:
            print("No results to save")
            return
        
        # Create summary data
        summary_data = []
        for result in results:
            summary_data.append({
                'filename': result['filename'],
                'fluid': result['fluid'],
                'is_mixture': result['is_mixture'],
                'x1': result['x1'],
                'x2': result['x2'],
                'TT101_mean/K': result['TT101_mean/K'],
                'n_points': result['n_points'],
                'polynomial_degree': result['polynomial_degree'],
                'r2_score': result['r2_score'],
                'mean_rel_abs_err_perc': result['mean_rel_abs_err_perc'],
                'std_rel_abs_err_perc': result['std_rel_abs_err_perc'],
                'rms_rel_abs_err_perc': result['rms_rel_abs_err_perc'],
                'mean_unc_perc': result['mean_unc_perc']
            })
        
        # Save summary as CSV
        summary_df = pd.DataFrame(summary_data)
        summary_output = os.path.join(self.output_dir, 'jt_coefficients_summary.csv')
        summary_df.to_csv(summary_output, index=False)
        
        # Save detailed results as JSON
        detailed_output = os.path.join(self.output_dir, 'jt_coefficients_detailed.json')
        
        # Convert numpy arrays to lists for JSON serialization
        json_results = []
        for result in results:
            json_result = result.copy()
            for key, value in json_result.items():
                if isinstance(value, np.ndarray):
                    json_result[key] = value.tolist()
                elif isinstance(value, (np.integer, np.floating)):
                    json_result[key] = float(value)
            json_results.append(json_result)
        
        with open(detailed_output, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        # Save individual CSV files for each measurement
        for result in results:
            csv_data = {
                'p/MPa': result['p/MPa'],
                'p_UNC/MPa': result['p_UNC/MPa'],
                'T/K': result['T/K'],
                'T_UNC/K': result['T_UNC/K'],
                'x1': result['x1_per_point'] if result['x1_per_point'] is not None else [np.nan] * result['n_points'],
                'JT_meas/(K/MPa)': result['JT_meas/(K/MPa)'],
                'JT_UNC/(K/MPa)': result['JT_UNC/(K/MPa)'],
                'JT_eos/(K/MPa)': result['JT_eos/(K/MPa)'],
                'mean_unc_perc': result['mean_unc_perc'],
                'rel_abs_err_perc': result['rel_abs_err_perc'],
                'used': result['used']
            }
            
            csv_df = pd.DataFrame(csv_data)
            csv_filename = result['filename'].replace('.csv', '_JT.csv')
            csv_output = os.path.join(self.output_dir, csv_filename)
            csv_df.to_csv(csv_output, index=False)
        
        return summary_df


def main():
    calculator = JTCoefficientCalculator()
    results = calculator.process_all_files()
    calculator.save_results(results)


if __name__ == "__main__":
    main()
