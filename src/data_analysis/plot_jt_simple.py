#!/usr/bin/env python3
"""
Simple plotting script for JT coefficient measurements.
Generates publication-quality plots without REFPROP theoretical calculations.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import Chebyshev
import warnings
warnings.filterwarnings('ignore')

# Add src path for imports
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'data_analysis'))
from get_git_root import get_git_root


def set_mpl():
    """Set matplotlib parameters for publication-quality plots."""
    fontsize = 16
    mpl.rcParams.update({
        "font.family": "Times New Roman",
        "mathtext.fontset": "dejavuserif",
        "font.size": fontsize,
        "axes.labelsize": fontsize,
        "axes.titlesize": fontsize,
        "legend.fontsize": fontsize,
        "xtick.top": True,
        "xtick.bottom": True,
        "ytick.left": True,
        "ytick.right": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.pad": 5,
        "ytick.major.pad": 5
    })


class JTPlotter:
    """Simple class for plotting JT coefficient results."""
    
    def __init__(self):
        """Initialize plotter with paths and load data."""
        self.git_root = get_git_root(os.getcwd())
        self.data_dir = os.path.join(self.git_root, 'data', 'derived_data', 'p_T_pairs')
        self.results_dir = os.path.join(self.git_root, 'data', 'derived_data', 'JT_coeffs')
        self.metadata_file = os.path.join(self.git_root, 'data', 'metadata', 'p_T_pairs.json')
        self.output_dir = os.path.join(self.git_root, 'figures')
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Load metadata and results
        self.metadata = self._load_metadata()
        self.results = self._load_results()
        
        # Set matplotlib style
        set_mpl()
    
    def _load_metadata(self):
        """Load metadata from JSON file."""
        try:
            with open(self.metadata_file, 'r') as f:
                metadata = json.load(f)
            return {item['filename']: item for item in metadata}
        except FileNotFoundError:
            print(f"Warning: Metadata file not found at {self.metadata_file}")
            return {}
    
    def _load_results(self):
        """Load JT coefficient results."""
        try:
            results_file = os.path.join(self.results_dir, 'jt_coefficients_detailed.json')
            with open(results_file, 'r') as f:
                results = json.load(f)
            return {item['filename']: item for item in results}
        except FileNotFoundError:
            print(f"Warning: Results file not found")
            return {}
    
    def _get_fluid_data(self, fluid_name):
        """Get all data for a specific fluid."""
        fluid_data = []
        for filename, result in self.results.items():
            if result['fluid'] == fluid_name:
                fluid_data.append(result)
        return fluid_data
    
    def plot_pure_fluids_with_errorbars(self, save_figure=True):
        """Plot pure fluid measurements with error bars."""
        # Get data for pure fluids
        nitrogen_data = self._get_fluid_data('Nitrogen')
        argon_data = self._get_fluid_data('Argon')
        helium_data = self._get_fluid_data('Helium')
        
        # Create subplots
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # Colors for different measurements
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
        
        # Plot Nitrogen
        self._plot_fluid_measurements(axes[0], nitrogen_data, 'Nitrogen', colors)
        axes[0].set_title('Nitrogen')
        axes[0].set_xlim(0, 12)
        axes[0].set_ylim(130, 170)
        
        # Plot Argon
        self._plot_fluid_measurements(axes[1], argon_data, 'Argon', colors)
        axes[1].set_title('Argon')
        axes[1].set_xlim(0, 12)
        axes[1].set_ylim(165, 185)
        
        # Plot Helium
        self._plot_fluid_measurements(axes[2], helium_data, 'Helium', colors)
        axes[2].set_title('Helium')
        axes[2].set_xlim(0, 12)
        axes[2].set_ylim(60, 145)
        
        # Set common labels
        for ax in axes:
            ax.set_xlabel('p / MPa')
            ax.set_ylabel('T / K')
            ax.grid(True, alpha=0.3)
        
        # Add legend to first subplot
        axes[0].legend(loc='upper right')
        
        plt.tight_layout()
        
        if save_figure:
            plt.savefig(os.path.join(self.output_dir, 'measurements_pures_errorbars.pdf'), 
                       bbox_inches='tight', dpi=300)
            plt.savefig(os.path.join(self.output_dir, 'measurements_pures_errorbars.png'), 
                       bbox_inches='tight', dpi=300)
        
        # plt.show()  # Disabled for non-interactive backend
        return fig, axes
    
    def plot_helium_neon_with_errorbars(self, save_figure=True):
        """Plot Helium-Neon mixture measurements with error bars."""
        # Get Helium-Neon data
        hene_data = self._get_fluid_data('Helium-Neon')
        
        # Group by composition
        composition_groups = {}
        for data in hene_data:
            comp = data.get('composition', 0)
            comp_key = f"{comp:.3f}"
            if comp_key not in composition_groups:
                composition_groups[comp_key] = []
            composition_groups[comp_key].append(data)
        
        # Sort by composition
        sorted_compositions = sorted(composition_groups.keys(), key=float)
        
        # Create figure with subplots for different compositions
        n_compositions = len(sorted_compositions)
        fig, axes = plt.subplots(2, (n_compositions + 1) // 2, figsize=(15, 10))
        if n_compositions <= 2:
            axes = axes.reshape(1, -1) if n_compositions == 2 else [axes]
        axes = axes.flatten()
        
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive']
        
        for i, comp_key in enumerate(sorted_compositions):
            if i >= len(axes):
                break
                
            comp_data = composition_groups[comp_key]
            self._plot_fluid_measurements(axes[i], comp_data, 'Helium-Neon', colors)
            
            # Set title with composition
            comp_value = float(comp_key)
            axes[i].set_title(f'He-Ne, x_He = {comp_value:.3f}')
            axes[i].set_xlim(0, 10)
            axes[i].set_ylim(50, 85)
            axes[i].set_xlabel('p / MPa')
            axes[i].set_ylabel('T / K')
            axes[i].grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(len(sorted_compositions), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        
        if save_figure:
            plt.savefig(os.path.join(self.output_dir, 'measurements_HeNe_errorbars.pdf'), 
                       bbox_inches='tight', dpi=300)
            plt.savefig(os.path.join(self.output_dir, 'measurements_HeNe_errorbars.png'), 
                       bbox_inches='tight', dpi=300)
        
        # plt.show()  # Disabled for non-interactive backend
        return fig, axes
    
    def _plot_fluid_measurements(self, ax, fluid_data, fluid_name, colors):
        """Plot measurements for a single fluid with error bars."""
        legend_added = False
        
        for i, data in enumerate(fluid_data):
            color = colors[i % len(colors)]
            
            # Get measurement data
            p_data = np.array(data['pressures_MPa'])
            T_data = np.array(data['temperatures_K'])
            p_unc = np.array(data['pressure_uncertainties_MPa'])
            T_unc = np.array(data['temperature_uncertainties_K'])
            
            # Plot measurement points with error bars
            ax.errorbar(p_data, T_data, xerr=p_unc, yerr=T_unc,
                       fmt='o', color=color, markersize=6, capsize=3,
                       label='measurements' if not legend_added else None)
            
            # Fit polynomial to data and plot
            if len(p_data) > 2:
                try:
                    poly_fit = Chebyshev.fit(p_data, T_data, min(len(p_data)-1, 4))
                    poly = Chebyshev(poly_fit.convert().coef)
                    
                    p_fit = np.linspace(p_data.min(), p_data.max(), 100)
                    T_fit = poly(p_fit)
                    
                    ax.plot(p_fit, T_fit, '--', color=color, linewidth=1.5,
                           label='fit' if not legend_added else None)
                except Exception as e:
                    print(f"Error fitting polynomial: {e}")
            
            legend_added = True
    
    def plot_jt_coefficients_comparison(self, save_figure=True):
        """Plot JT coefficients with error bars comparing measured vs theoretical."""
        # Get all fluid data
        all_fluids = ['Nitrogen', 'Argon', 'Helium', 'Helium-Neon', 'Helium-Nitrogen']
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        for i, fluid_name in enumerate(all_fluids):
            if i >= len(axes):
                break
                
            fluid_data = self._get_fluid_data(fluid_name)
            self._plot_jt_coefficients(axes[i], fluid_data, fluid_name)
            
            axes[i].set_title(fluid_name)
            axes[i].set_xlabel('p / MPa')
            axes[i].set_ylabel('μ_JT / (K/MPa)')
            axes[i].grid(True, alpha=0.3)
        
        # Hide unused subplot
        if len(all_fluids) < len(axes):
            axes[-1].set_visible(False)
        
        plt.tight_layout()
        
        if save_figure:
            plt.savefig(os.path.join(self.output_dir, 'jt_coefficients_comparison.pdf'), 
                       bbox_inches='tight', dpi=300)
            plt.savefig(os.path.join(self.output_dir, 'jt_coefficients_comparison.png'), 
                       bbox_inches='tight', dpi=300)
        
        # plt.show()  # Disabled for non-interactive backend
        return fig, axes
    
    def _plot_jt_coefficients(self, ax, fluid_data, fluid_name):
        """Plot JT coefficients for a single fluid."""
        colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
        
        for i, data in enumerate(fluid_data):
            color = colors[i % len(colors)]
            
            p_data = np.array(data['pressures_MPa'])
            jt_measured = np.array(data['jt_measured_K_per_MPa'])
            jt_theoretical = np.array(data['jt_theoretical_K_per_MPa'])
            jt_uncertainty = np.array(data['jt_uncertainty_K_per_MPa'])
            
            # Plot measured JT coefficients with error bars
            ax.errorbar(p_data, jt_measured, yerr=jt_uncertainty,
                       fmt='o', color=color, markersize=6, capsize=3,
                       label=f'Measured ({i+1})' if len(fluid_data) > 1 else 'Measured')
            
            # Plot theoretical JT coefficients
            ax.plot(p_data, jt_theoretical, 's', color=color, markersize=4,
                   label=f'EOS ({i+1})' if len(fluid_data) > 1 else 'EOS')
        
        ax.legend()


def main():
    """Main function to generate all plots."""
    print("Generating JT coefficient plots...")
    
    # Create plotter instance
    plotter = JTPlotter()
    
    # Generate plots
    print("1. Plotting pure fluids with error bars...")
    plotter.plot_pure_fluids_with_errorbars()
    
    print("2. Plotting Helium-Neon mixtures with error bars...")
    plotter.plot_helium_neon_with_errorbars()
    
    print("3. Plotting JT coefficients comparison...")
    plotter.plot_jt_coefficients_comparison()
    
    print(f"Plots saved to: {plotter.output_dir}")


if __name__ == "__main__":
    main()