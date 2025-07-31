#!/usr/bin/env python3
"""
Extract file information from p_T_pairs directory.
Creates a JSON file with filename, fluid name, and mean TT101/K for each file.
"""

import os
import json
import pandas as pd
from get_git_root import get_git_root


def parse_fluid_from_filename(filename):
    """Parse fluid name from filename."""
    # Remove .csv extension
    name = filename.replace('.csv', '')
    
    if 'Helium-Neon' in name:
        return 'Helium-Neon'
    elif 'Helium-Nitrogen' in name:
        return 'Helium-Nitrogen'
    elif 'Nitrogen' in name:
        return 'Nitrogen'
    elif 'Argon' in name:
        return 'Argon'
    elif 'Helium' in name:
        return 'Helium'
    else:
        return 'Unknown'


def main():
    """Extract file information and create JSON file."""
    git_root = get_git_root(os.getcwd())
    data_dir = os.path.join(git_root, 'data', 'derived_data', 'p_T_pairs')
    
    # Get all CSV files
    csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
    csv_files.sort()
    
    file_info = []
    
    print(f"Processing {len(csv_files)} files...")
    
    for csv_file in csv_files:
        filepath = os.path.join(data_dir, csv_file)
        
        try:
            # Read the CSV file
            df = pd.read_csv(filepath)
            
            # Extract information
            fluid_name = parse_fluid_from_filename(csv_file)
            mean_tt101 = df['TT101/K'].mean()
            
            # Add to list
            file_info.append({
                'filename': csv_file,
                'fluid': fluid_name,
                'mean_TT101_K': mean_tt101
            })
            
            print(f"  {csv_file}: {fluid_name}, mean TT101/K = {mean_tt101:.3f}")
            
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")
    
    # Save to JSON file
    output_file = os.path.join(data_dir, 'file_info.json')
    with open(output_file, 'w') as f:
        json.dump(file_info, f, indent=2)
    
    print(f"\nFile information saved to: {output_file}")
    print(f"Total files processed: {len(file_info)}")


if __name__ == "__main__":
    main()