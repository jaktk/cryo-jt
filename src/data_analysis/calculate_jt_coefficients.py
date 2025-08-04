import os
import json
import numpy as np
import pandas as pd
from scipy import optimize
from numpy.polynomial.chebyshev import Chebyshev
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings('ignore')

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'data_analysis'))
from FluidProps import FluidProps
from get_git_root import get_git_root


class JTCoefficientCalculator:
    """ Calculate Joule-Thomson coefficients from p_T_pairs data """
    
    def __init__(self):
        """ Initialize the calculator with paths and create output directory """
        self.git_root = get_git_root(os.getcwd())
        self.data_dir = os.path.join(self.git_root, 'data', 'derived_data', 'p_T_pairs')
        self.metadata_file = os.path.join(self.git_root, 'data', 'metadata', 'p_T_pairs.json')
        self.output_dir = os.path.join(self.git_root, 'data', 'derived_data', 'JT_coeffs')
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load metadata
        self.metadata = self._load_metadata()
        
    def _load_metadata(self):
        """Load metadata from JSON file."""
        try:
            with open(self.metadata_file, 'r') as f:
                metadata = json.load(f)
            
            # Convert to dictionary keyed by filename for easy lookup
            metadata_dict = {}
            for item in metadata:
                metadata_dict[item['filename']] = item
            
            return metadata_dict
        except FileNotFoundError:
            print(f"Warning: Metadata file not found at {self.metadata_file}")
            return {}
    
    def _parse_fluid_info(self, filename):
        """ Parse fluid type and composition from metadata """
        # Get metadata for this file
        metadata = self.metadata.get(filename, {})
        
        if not metadata:
            print(f"Warning: No metadata found for {filename}")
            return None, False
        
        fluid_name = metadata.get('fluid', '')
        
        # Determine fluid type and mixture status from metadata
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
        
                # Stop if R^2 is sufficiently high
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
            # Read data
            df = pd.read_csv(filepath)
            
            # Get metadata
            metadata = self.metadata.get(filename, {})
            
            # Parse fluid information from metadata
            fluid_type, is_mixture = self._parse_fluid_info(filename)
            
            if fluid_type is None:
                print(f"  Skipping {filename} - could not parse fluid info")
                return None
            
            # Extract pressure and temperature data
            pressures = df['PT102/MPa'].values # Downstream pressure
            temperatures = df['TT102/K'].values # Downstream temperature
            
            # Get measurement uncertainties
            p_uncertainties = df['PT102/MPa_EXP_UNC'].values
            t_uncertainties = df['TT102/K_EXP_UNC'].values
            
            # Skip if insufficient data
            if len(pressures) < 3:
                print(f"  Insufficient data points: {len(pressures)}")
                return None
            
            # Initialize FluidProps
            fluid_props = FluidProps(fluid_type)
            
            # For mixtures, get composition
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
                'composition': composition,
                'n_points': len(pressures),
                'polynomial_degree': degree,
                'r2_score': r2,
                'mean_TT101_K': metadata.get('mean_TT101_K', np.nan),
                'pressures_MPa': pressures,
                'temperatures_K': temperatures,
                'pressure_uncertainties_MPa': p_uncertainties,
                'temperature_uncertainties_K': t_uncertainties,
                'jt_measured_K_per_MPa': jt_measured,
                'jt_theoretical_K_per_MPa': jt_theoretical,
                'jt_uncertainty_K_per_MPa': jt_uncertainty,
                'relative_abs_error_percent': relative_abs_error,
                'mean_relative_abs_error_percent': np.nanmean(relative_abs_error),
                'std_relative_abs_error_percent': np.nanstd(relative_abs_error),
                'rms_relative_abs_error_percent': np.sqrt(np.nanmean(relative_abs_error**2)),
                'mean_uncertainty_percent': np.nanmean(jt_uncertainty / np.abs(jt_measured) * 100)
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
                'composition': result['composition'],
                'mean_TT101_K': result['mean_TT101_K'],
                'n_points': result['n_points'],
                'polynomial_degree': result['polynomial_degree'],
                'r2_score': result['r2_score'],
                'mean_relative_abs_error_percent': result['mean_relative_abs_error_percent'],
                'std_relative_abs_error_percent': result['std_relative_abs_error_percent'],
                'rms_relative_abs_error_percent': result['rms_relative_abs_error_percent'],
                'mean_uncertainty_percent': result['mean_uncertainty_percent']
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
                'pressure_MPa': result['pressures_MPa'],
                'temperature_K': result['temperatures_K'],
                'pressure_uncertainty_MPa': result['pressure_uncertainties_MPa'],
                'temperature_uncertainty_K': result['temperature_uncertainties_K'],
                'jt_measured_K_per_MPa': result['jt_measured_K_per_MPa'],
                'jt_theoretical_K_per_MPa': result['jt_theoretical_K_per_MPa'],
                'jt_uncertainty_K_per_MPa': result['jt_uncertainty_K_per_MPa'],
                'relative_abs_error_percent': result['relative_abs_error_percent']
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
