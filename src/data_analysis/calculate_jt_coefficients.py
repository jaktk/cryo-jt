#!/usr/bin/env python3
"""
Calculate Joule-Thomson coefficients from p_T_pairs data.

This script processes all files in data/derived_data/p_T_pairs/ and calculates
JT coefficients by fitting Chebyshev polynomials to the isenthalpic data and 
taking derivatives. Uses metadata from data/metadata/p_T_pairs.json.

Implements full uncertainty propagation through polynomial fitting and differentiation.

Author: Based on PhD thesis Chapter 4 methodology
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import optimize
from numpy.polynomial.chebyshev import Chebyshev
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings('ignore')

# Add src path for imports
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'data_analysis'))
from FluidProps import FluidProps
from get_git_root import get_git_root


class JTCoefficientCalculator:
    """Calculate Joule-Thomson coefficients from p_T_pairs data."""
    
    def __init__(self):
        """Initialize the calculator with paths and create output directory."""
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
        """Parse fluid type and composition from metadata."""
        # Get metadata for this file
        metadata = self.metadata.get(filename, {})
        
        if not metadata:
            print(f"Warning: No metadata found for {filename}")
            return None, False
        
        fluid_name = metadata.get('fluid', '')
        
        # Determine fluid type and mixture status from metadata
        if fluid_name == 'Helium-Neon':
            fluid_type = ('Helium', 'Neon')
            mixture = True
        elif fluid_name == 'Nitrogen-Helium':
            fluid_type = ('Nitrogen', 'Helium')
            mixture = True
        elif fluid_name == 'Nitrogen':
            fluid_type = 'Nitrogen'
            mixture = False
        elif fluid_name == 'Argon':
            fluid_type = 'Argon'
            mixture = False
        elif fluid_name == 'Helium':
            fluid_type = 'Helium'
            mixture = False
        else:
            print(f"Warning: Unknown fluid type '{fluid_name}' in metadata for {filename}")
            return None, False
        
        return fluid_type, mixture
    
    def _fit_polynomial_weighted(self, pressures, temperatures, p_uncertainties, t_uncertainties, max_degree=6):
        """
        Fit optimal Chebyshev polynomial to p-T data with weighted least squares.
        
        Uses measurement uncertainties to weight the fitting process.
        """
        n_points = len(pressures)
        max_degree = min(max_degree, n_points - 1)
        
        # Calculate weights from uncertainties
        # Weight = 1/σ² for each measurement
        weights = 1.0 / (t_uncertainties**2 + (p_uncertainties * np.gradient(temperatures, pressures))**2)
        weights = weights / np.sum(weights)  # Normalize weights
        
        best_degree = 1
        best_r2 = -np.inf
        best_fit = None
        best_cov = None
        
        for degree in range(1, max_degree + 1):
            try:
                # Fit weighted Chebyshev polynomial
                poly_fit = Chebyshev.fit(pressures, temperatures, degree, w=weights)
                poly = Chebyshev(poly_fit.convert().coef)
                
                # Calculate R²
                fitted_temps = poly(pressures)
                r2 = r2_score(temperatures, fitted_temps)
                
                # Store best fit
                if r2 > best_r2:
                    best_r2 = r2
                    best_degree = degree
                    best_fit = poly
                    
                    # Calculate covariance matrix for uncertainty propagation
                    best_cov = self._calculate_covariance_matrix(
                        pressures, temperatures, poly, weights, degree
                    )
                
                # Stop if R² is sufficiently high
                if r2 > 0.999:
                    break
                    
            except Exception as e:
                print(f"Error fitting degree {degree}: {e}")
                continue
        
        return best_fit, best_degree, best_r2, best_cov
    
    def _calculate_covariance_matrix(self, pressures, temperatures, polynomial, weights, degree):
        """
        Calculate covariance matrix for polynomial coefficients.
        
        This is used for uncertainty propagation through the polynomial derivative.
        """
        try:
            # Create design matrix for Chebyshev polynomials
            n_points = len(pressures)
            n_coeffs = degree + 1
            
            # Design matrix A where A[i,j] = T_j(p_i) (Chebyshev polynomial)
            A = np.zeros((n_points, n_coeffs))
            
            # Map pressure to [-1, 1] domain for Chebyshev polynomials
            p_min, p_max = np.min(pressures), np.max(pressures)
            p_scaled = 2 * (pressures - p_min) / (p_max - p_min) - 1
            
            for i in range(n_points):
                for j in range(n_coeffs):
                    if j == 0:
                        A[i, j] = 1.0
                    elif j == 1:
                        A[i, j] = p_scaled[i]
                    else:
                        # Chebyshev recurrence relation
                        A[i, j] = 2 * p_scaled[i] * A[i, j-1] - A[i, j-2]
            
            # Weight the design matrix
            W = np.diag(weights)
            
            # Calculate covariance matrix: (A^T W A)^-1
            AtWA = A.T @ W @ A
            
            # Add small regularization to prevent singularity
            reg_factor = 1e-12 * np.trace(AtWA) / n_coeffs
            AtWA += reg_factor * np.eye(n_coeffs)
            
            cov_matrix = np.linalg.inv(AtWA)
            
            return cov_matrix
            
        except Exception as e:
            print(f"Error calculating covariance matrix: {e}")
            # Return identity matrix as fallback
            return np.eye(degree + 1)
    
    def _calculate_jt_coefficients_with_uncertainty(self, pressures, polynomial, covariance_matrix):
        """
        Calculate JT coefficients with full uncertainty propagation.
        
        Uses the covariance matrix of polynomial coefficients to propagate
        uncertainties through the derivative calculation.
        """
        # Get derivative of polynomial (this gives dT/dp)
        poly_deriv = polynomial.deriv(m=1)
        
        # Calculate JT coefficients at measurement points
        jt_coefficients = poly_deriv(pressures)
        
        # Calculate uncertainties using covariance matrix
        jt_uncertainties = np.zeros_like(jt_coefficients)
        
        try:
            # Get derivative coefficients
            deriv_coeffs = poly_deriv.coef
            
            # Map pressure to [-1, 1] domain for Chebyshev polynomials
            p_min, p_max = np.min(pressures), np.max(pressures)
            p_scaled = 2 * (pressures - p_min) / (p_max - p_min) - 1
            
            for i, p_val in enumerate(p_scaled):
                # Calculate derivative basis functions at this pressure
                deriv_basis = np.zeros(len(deriv_coeffs))
                
                if len(deriv_coeffs) > 0:
                    deriv_basis[0] = 1.0
                if len(deriv_coeffs) > 1:
                    deriv_basis[1] = p_val
                
                for j in range(2, len(deriv_coeffs)):
                    deriv_basis[j] = 2 * p_val * deriv_basis[j-1] - deriv_basis[j-2]
                
                # Scale factor for domain transformation
                scale_factor = 2 / (p_max - p_min)
                deriv_basis *= scale_factor
                
                # Propagate uncertainty: σ² = basis^T * Cov * basis
                # Only use the part of covariance matrix corresponding to derivative coefficients
                n_deriv = len(deriv_coeffs)
                if n_deriv <= covariance_matrix.shape[0]:
                    cov_sub = covariance_matrix[:n_deriv, :n_deriv]
                    variance = deriv_basis.T @ cov_sub @ deriv_basis
                    jt_uncertainties[i] = np.sqrt(max(0, variance))
                else:
                    # Fallback: use simplified uncertainty estimation
                    jt_uncertainties[i] = np.abs(jt_coefficients[i]) * 0.05  # 5% uncertainty
                    
        except Exception as e:
            print(f"Error in uncertainty propagation: {e}")
            # Fallback to simplified uncertainty
            jt_uncertainties = np.abs(jt_coefficients) * 0.05  # 5% uncertainty
        
        return jt_coefficients, jt_uncertainties
    
    def _calculate_theoretical_jt(self, pressures, temperatures, fluid_props, composition=None):
        """Calculate theoretical JT coefficients using FluidProps."""
        jt_theoretical = np.zeros_like(pressures)
        
        for i, (p, t) in enumerate(zip(pressures, temperatures)):
            try:
                if composition is not None:
                    # For mixtures, set composition
                    fluid_props.set_composition_from_1st_fraction(composition)
                
                jt_theoretical[i] = fluid_props.get_JT_coefficient(p, t)
                
            except Exception as e:
                print(f"Error calculating theoretical JT for p={p}, T={t}: {e}")
                jt_theoretical[i] = np.nan
        
        return jt_theoretical
    
    def _monte_carlo_uncertainty(self, pressures, temperatures, p_uncertainties, t_uncertainties, 
                                polynomial, n_samples=1000):
        """
        Monte Carlo uncertainty estimation for JT coefficients.
        
        This provides an independent check of the analytical uncertainty propagation.
        """
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
        """Process a single p_T_pairs file to calculate JT coefficients."""
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
            pressures = df['PT102/MPa'].values  # Downstream pressure
            temperatures = df['TT102/K'].values  # Downstream temperature
            
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
            poly, degree, r2, cov_matrix = self._fit_polynomial_weighted(
                pressures, temperatures, p_uncertainties, t_uncertainties
            )
            
            if poly is None:
                print(f"  Failed to fit polynomial")
                return None
            
            # Calculate JT coefficients with full uncertainty propagation
            jt_measured, jt_uncertainty_analytical = self._calculate_jt_coefficients_with_uncertainty(
                pressures, poly, cov_matrix
            )
            
            # Calculate Monte Carlo uncertainty for comparison
            jt_uncertainty_mc = self._monte_carlo_uncertainty(
                pressures, temperatures, p_uncertainties, t_uncertainties, poly
            )
            
            # Use the larger of analytical and Monte Carlo uncertainties
            jt_uncertainty = np.maximum(jt_uncertainty_analytical, jt_uncertainty_mc)
            
            # Calculate theoretical JT coefficients
            jt_theoretical = self._calculate_theoretical_jt(pressures, temperatures, fluid_props, composition)
            
            # Calculate relative errors
            relative_error = np.where(jt_theoretical != 0, 
                                    (jt_measured - jt_theoretical) / jt_theoretical * 100,
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
                'jt_uncertainty_analytical_K_per_MPa': jt_uncertainty_analytical,
                'jt_uncertainty_mc_K_per_MPa': jt_uncertainty_mc,
                'jt_uncertainty_K_per_MPa': jt_uncertainty,
                'relative_error_percent': relative_error,
                'mean_relative_error_percent': np.nanmean(relative_error),
                'std_relative_error_percent': np.nanstd(relative_error),
                'rms_relative_error_percent': np.sqrt(np.nanmean(relative_error**2)),
                'mean_uncertainty_percent': np.nanmean(jt_uncertainty / np.abs(jt_measured) * 100)
            }
            
            print(f"  Success: {len(pressures)} points, degree={degree}, R²={r2:.6f}")
            print(f"  Mean relative error: {np.nanmean(relative_error):.2f}%")
            print(f"  Mean uncertainty: {np.nanmean(jt_uncertainty / np.abs(jt_measured) * 100):.2f}%")
            
            return results
            
        except Exception as e:
            print(f"  Error processing {filename}: {e}")
            return None
    
    def process_all_files(self):
        """Process all files in the p_T_pairs directory."""
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
        """Save results to files."""
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
                'mean_relative_error_percent': result['mean_relative_error_percent'],
                'std_relative_error_percent': result['std_relative_error_percent'],
                'rms_relative_error_percent': result['rms_relative_error_percent'],
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
                'jt_uncertainty_analytical_K_per_MPa': result['jt_uncertainty_analytical_K_per_MPa'],
                'jt_uncertainty_mc_K_per_MPa': result['jt_uncertainty_mc_K_per_MPa'],
                'jt_uncertainty_K_per_MPa': result['jt_uncertainty_K_per_MPa'],
                'relative_error_percent': result['relative_error_percent']
            }
            
            csv_df = pd.DataFrame(csv_data)
            csv_filename = result['filename'].replace('.csv', '_JT.csv')
            csv_output = os.path.join(self.output_dir, csv_filename)
            csv_df.to_csv(csv_output, index=False)
        
        print(f"\\nResults saved to: {self.output_dir}")
        print(f"- Summary table: jt_coefficients_summary.csv")
        print(f"- Detailed results: jt_coefficients_detailed.json")
        print(f"- Individual CSV files: *_JT.csv")
        
        return summary_df
    
    def print_summary(self, summary_df):
        """Print analysis summary."""
        print(f"\\n{'='*60}")
        print("JOULE-THOMSON COEFFICIENT ANALYSIS SUMMARY")
        print(f"{'='*60}")
        
        print(f"Total files processed: {len(summary_df)}")
        
        # Group by fluid type
        fluid_counts = summary_df.groupby('fluid').size()
        print(f"\\nFluid distribution:")
        for fluid, count in fluid_counts.items():
            print(f"  {fluid}: {count} files")
        
        # Temperature and pressure ranges
        print(f"\\nExperimental conditions:")
        print(f"  Temperature range: {summary_df['mean_TT101_K'].min():.1f}K to {summary_df['mean_TT101_K'].max():.1f}K")
        
        # Data quality metrics
        print(f"\\nData quality metrics:")
        print(f"  Mean polynomial degree: {summary_df['polynomial_degree'].mean():.1f}")
        print(f"  Mean R² score: {summary_df['r2_score'].mean():.6f}")
        print(f"  Mean relative error: {summary_df['mean_relative_error_percent'].mean():.2f}%")
        print(f"  RMS relative error: {summary_df['rms_relative_error_percent'].mean():.2f}%")
        print(f"  Mean uncertainty: {summary_df['mean_uncertainty_percent'].mean():.2f}%")
        
        # Best and worst cases
        best_idx = summary_df['rms_relative_error_percent'].idxmin()
        worst_idx = summary_df['rms_relative_error_percent'].idxmax()
        
        print(f"\\nBest measurement:")
        print(f"  File: {summary_df.loc[best_idx, 'filename']}")
        print(f"  Fluid: {summary_df.loc[best_idx, 'fluid']}")
        print(f"  RMS error: {summary_df.loc[best_idx, 'rms_relative_error_percent']:.2f}%")
        
        print(f"\\nWorst measurement:")
        print(f"  File: {summary_df.loc[worst_idx, 'filename']}")
        print(f"  Fluid: {summary_df.loc[worst_idx, 'fluid']}")
        print(f"  RMS error: {summary_df.loc[worst_idx, 'rms_relative_error_percent']:.2f}%")
        
        print(f"\\n{'='*60}")


def main():
    """Main function to run the JT coefficient calculation."""
    print("Joule-Thomson Coefficient Calculator")
    print("Based on PhD thesis Chapter 4 methodology")
    print("Full uncertainty propagation implementation")
    print("="*50)
    
    # Initialize calculator
    calculator = JTCoefficientCalculator()
    
    # Process all files
    results = calculator.process_all_files()
    
    if not results:
        print("No files processed successfully")
        return
    
    # Save results
    summary_df = calculator.save_results(results)
    
    # Print summary
    calculator.print_summary(summary_df)


if __name__ == "__main__":
    main()