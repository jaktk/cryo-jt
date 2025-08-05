import os
import json
import numpy as np
import pandas as pd
from scipy import optimize
from numpy.polynomial.chebyshev import Chebyshev
from sklearn.metrics import r2_score
from FluidProps import FluidProps
from get_git_root import get_git_root


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

    def get_polynomial_fit_uncert(self, T):
        """ Uncertainty in mK from fitting Chebyshev polynomial """
        d = {
        'X93303': [{'Tmin': 1.4 , 'Tmax': 14.1 , 'N': 31, 'n': 9, 'DTrms': 0.93},
                  {'Tmin': 14.1, 'Tmax': 80.0 , 'N': 31, 'n': 7, 'DTrms': 1.60},
                  {'Tmin': 80.0, 'Tmax': 325.0, 'N': 32, 'n': 8, 'DTrms': 5.84}],
        'X115143': [{'Tmin': 20.0, 'Tmax': 95.0, 'N': 28, 'n': 7, 'DTrms': 0.91},
                   {'Tmin': 95.0, 'Tmax': 325.0, 'N': 29, 'n': 9, 'DTrms': 3.45}],
        'X115888': [{'Tmin': 20.0, 'Tmax': 95.3, 'N': 28, 'n': 7, 'DTrms': 1.17},
                   {'Tmin': 95.3, 'Tmax': 325.0, 'N': 29, 'n': 8, 'DTrms': 4.06}]
        }
        for _ in d[self.serNum]:
            if T >= _['Tmin'] and T < _['Tmax']:
                sigm2 = _['N'] / (_['N'] - _['n']) * _['DTrms']**2
                break
            else:
                sigm2 = 0
        return 2 * sigm2**0.5


class TempUncertainty:
    sensor = ''
    def __init__(self, sensor):
        assert sensor in ['X93303','X115143','X115888']
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
        'X115888': [{'Tmin': 20.0, 'Tmax': 95.3, 'N': 28, 'n': 7, 'DTrms': 1.17},
                   {'Tmin': 95.3, 'Tmax': 325.0, 'N': 29, 'n': 8, 'DTrms': 4.06}]
        }

        for i in d[self.sensor]:
            if T >= i['Tmin'] and T < i['Tmax']:
                sigm2 = i['N'] / (i['N'] - i['n']) * i['DTrms']**2
                break
            else:
                sigm2 = 0
        return 2 * sigm2**0.5

    def __call__(self, T):
        """ return combined temperature uncertainty in K """
        cabtr = self.get_cabtr_uncertainty(T)
        cernox = self.get_cernox_uncertainty(T)
        poly = self.get_polynomial_fit_uncert(T)
        return (cabtr + cernox + poly) * 1e-3


class JTCoefficientCalculator:
    def __init__(self):
        self.git_root = get_git_root(os.getcwd())
        self.data_dir = os.path.join(self.git_root, 'data', 'derived_data', 'p_T_pairs')
        self.metadata_file = os.path.join(self.git_root, 'data', 'metadata', 'p_T_pairs.json')
        self.output_dir = os.path.join(self.git_root, 'data', 'derived_data', 'JT_coeffs')
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
    
    def _fit_polynomial_weighted(self, pressures, temperatures, p_uncertainties, t_uncertainties, max_degree=6):
        """
        Fit optimal Chebyshev polynomial to p-T data with weighted least squares.
        Uses measurement uncertainties to weight the fitting process.
        """
        n_points = len(pressures)
        max_degree = min(max_degree, n_points - 1)
        
        # Calculate weights from uncertainties
        # Weight = 1/σ^2 for each measurement
        weights = 1.0 / (t_uncertainties**2 + (p_uncertainties * np.gradient(temperatures, pressures))**2)
        weights = weights / np.sum(weights)  # Normalize weights
        
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
    
    def _calculate_theoretical_jt(self, pressures, temperatures, fluid_props, composition=None):
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
    
    def _monte_carlo_uncertainty(self, pressures, temperatures, p_uncertainties, t_uncertainties, polynomial, n_samples=1000):
        """ Monte Carlo uncertainty estimation for JT coefficients """
        try:
            jt_samples = []
            
            for _ in range(n_samples):
                # Generate random samples within measurement uncertainties
                p_sample = np.random.normal(pressures, p_uncertainties)
                t_sample = np.random.normal(temperatures, t_uncertainties)
                
                # Fit polynomial to perturbed data
                try:
                    poly_fit = Chebyshev.fit(p_sample, t_sample, polynomial.degree())
                    poly_sample = Chebyshev(poly_fit.convert().coef)
                    
                    # Calculate JT coefficients
                    poly_deriv = poly_sample.deriv(m=1)
                    jt_sample = poly_deriv(pressures)
                    
                    jt_samples.append(jt_sample)
                    
                except:
                    continue
            
            if jt_samples:
                jt_samples = np.array(jt_samples)
                jt_mc_uncertainty = np.std(jt_samples, axis=0)
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
            
            # Get measurement uncertainties (max value of sensor uncertainty and standard expanded uncertainty from a set of measurements)
            p_uncertainties = np.array([max(p, 1.96 * 0.01/100 * 13.7) for p in df['PT102/MPa_EXP_UNC'].values])
            comb_T_unc = TempUncertainty("X93303")
            t_uncertainties = np.array([max(t, comb_T_unc(t)) for t in df['TT102/K_EXP_UNC'].values])
            
            # Skip if insufficient data
            if len(pressures) < 3:
                print(f"  Insufficient data points: {len(pressures)}")
                return None
            
            fluid_props = FluidProps(fluid_type)
            
            # Get composition, if mixture
            composition = None
            if is_mixture:
                if 'x1' in df.columns and not df['x1'].isna().all():
                    composition = df['x1'].mean()
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
            
            # Calculate relative errors
            relative_abs_error = np.where(jt_theoretical != 0, 
                                    np.abs((jt_measured - jt_theoretical) / jt_theoretical * 100),
                                    np.nan)
            
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
                'p_UNC/MPa': p_uncertainties,
                'T_UNC/K': t_uncertainties,
                'JT_meas/(K/MPa)': jt_measured,
                'JT_eos/(K/MPa)': jt_theoretical,
                'JT_UNC/(K/MPa)': jt_uncertainty,
                'rel_abs_err_perc': relative_abs_error,
                'mean_rel_abs_err_perc': np.nanmean(relative_abs_error),
                'std_rel_abs_err_perc': np.nanstd(relative_abs_error),
                'rms_rel_abs_err_perc': np.sqrt(np.nanmean(relative_abs_error**2)),
                'mean_unc_perc': np.nanmean(jt_uncertainty / np.abs(jt_measured) * 100)
            }
            
            print(f"  Success: {len(pressures)} points, degree={degree}, R²={r2:.6f}")
            print(f"  Mean relative error: {np.nanmean(relative_abs_error):.2f}%")
            print(f"  Mean uncertainty: {np.nanmean(jt_uncertainty / np.abs(jt_measured) * 100):.2f}%")
            
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
                'T/K': result['T/K'],
                'p_UNC/MPa': result['p_UNC/MPa'],
                'T_UNC/K': result['T_UNC/K'],
                'JT_meas/(K/MPa)': result['JT_meas/(K/MPa)'],
                'JT_eos/(K/MPa)': result['JT_eos/(K/MPa)'],
                'JT_UNC/(K/MPa)': result['JT_UNC/(K/MPa)'],
                'mean_rel_abs_err_perc': result['mean_rel_abs_err_perc']
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
