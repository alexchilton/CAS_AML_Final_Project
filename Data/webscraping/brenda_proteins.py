#!/usr/bin/env python3
"""
Script to extract pH data from BRENDA database tarball
Handles .tar.gz files from BRENDA
"""

import tarfile
import requests
import pandas as pd
import json
import time
import os
import re
from pathlib import Path

def extract_brenda_tarball(tarball_path="brenda_2024_1.txt.tar.gz"):
    """Extract BRENDA tarball"""
    print(f"Extracting {tarball_path}...")
    
    try:
        with tarfile.open(tarball_path, "r:gz") as tar:
            # Extract to current directory
            tar.extractall()
            # Get the extracted filename
            members = tar.getnames()
            txt_files = [m for m in members if m.endswith('.txt')]
            if txt_files:
                return txt_files[0]
    except Exception as e:
        print(f"Error extracting tarball: {e}")
    
    return None

def parse_brenda_flatfile(filename):
    """Parse BRENDA flatfile format for pH data"""
    if not os.path.exists(filename):
        print(f"File {filename} not found!")
        return None
    
    print(f"Parsing {filename}...")
    enzymes = {}
    current_ec = None
    
    with open(filename, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            
            # EC number line (starts with ID)
            if line.startswith("ID\t"):
                current_ec = line.split('\t')[1]
                enzymes[current_ec] = {
                    'name': '',
                    'ph_optimum': [],
                    'ph_range': [],
                    'ph_stability': []
                }
            
            # Protein name (starts with PR)
            elif line.startswith("PR\t") and current_ec:
                enzymes[current_ec]['name'] = line.split('\t', 1)[1]
            
            # pH optimum (starts with PH_OPT)
            elif line.startswith("PH_OPTIMUM\t") and current_ec:
                ph_data = parse_ph_entry(line)
                if ph_data:
                    enzymes[current_ec]['ph_optimum'].extend(ph_data)
            
            # pH range (starts with PH_RANGE)
            elif line.startswith("PH_RANGE\t") and current_ec:
                ph_data = parse_ph_entry(line)
                if ph_data:
                    enzymes[current_ec]['ph_range'].extend(ph_data)
            
            # pH stability (starts with PH_STABILITY)
            elif line.startswith("PH_STABILITY\t") and current_ec:
                ph_data = parse_ph_entry(line)
                if ph_data:
                    enzymes[current_ec]['ph_stability'].extend(ph_data)
    
    print(f"Parsed {len(enzymes)} enzymes from BRENDA")
    return enzymes

def parse_ph_entry(line):
    """Parse a pH data line from BRENDA"""
    try:
        # Remove the field name
        parts = line.split('\t', 1)
        if len(parts) < 2:
            return []
        
        content = parts[1]
        entries = []
        
        # Split multiple entries (separated by #)
        items = content.split('#')
        
        for item in items:
            if not item.strip():
                continue
            
            # Extract pH value using regex
            ph_matches = re.findall(r'(\d+(?:\.\d+)?)', item)
            
            if ph_matches:
                # Extract organism info if present
                organism = ""
                if '{' in item and '}' in item:
                    organism_match = re.search(r'\{([^}]+)\}', item)
                    if organism_match:
                        organism = organism_match.group(1)
                
                # Extract reference if present
                reference = ""
                if '<' in item and '>' in item:
                    ref_match = re.search(r'<([^>]+)>', item)
                    if ref_match:
                        reference = ref_match.group(1)
                
                # For pH ranges, take the midpoint
                if len(ph_matches) == 2 and '-' in item:
                    ph_value = (float(ph_matches[0]) + float(ph_matches[1])) / 2
                else:
                    ph_value = float(ph_matches[0])
                
                entries.append({
                    'value': ph_value,
                    'organism': organism,
                    'reference': reference,
                    'original': item.strip()
                })
    
    except Exception as e:
        print(f"Error parsing pH entry: {e}")
    
    return entries

def match_brenda_to_pdb(brenda_data, pdb_structures):
    """Match BRENDA pH data with PDB structures"""
    matched_data = []
    
    for pdb_id, pdb_info in pdb_structures.items():
        ec_number = pdb_info.get('ec_number')
        
        if ec_number and ec_number in brenda_data:
            brenda_info = brenda_data[ec_number]
            
            # Get pH values in order of preference
            ph_entries = []
            
            # 1. First try pH optimum
            if brenda_info['ph_optimum']:
                ph_entries = brenda_info['ph_optimum']
                source = 'pH_optimum'
            # 2. Then try pH range (less accurate)
            elif brenda_info['ph_range']:
                ph_entries = brenda_info['ph_range']
                source = 'pH_range'
            # 3. Finally try pH stability
            elif brenda_info['ph_stability']:
                ph_entries = brenda_info['ph_stability']
                source = 'pH_stability'
            
            if ph_entries:
                # Take the first value or median if multiple
                ph_values = [entry['value'] for entry in ph_entries if 'value' in entry]
                
                if ph_values:
                    median_ph = sorted(ph_values)[len(ph_values)//2]
                    
                    matched_data.append({
                        'pdb_id': pdb_id,
                        'ec_number': ec_number,
                        'enzyme_name': pdb_info.get('enzyme_name', ''),
                        'ph_optimum': median_ph,
                        'ph_source': f'BRENDA_{source}',
                        'brenda_name': brenda_info['name'],
                        'ph_values_count': len(ph_values),
                        'organisms': ', '.join(set(e['organism'] for e in ph_entries if e['organism'])),
                        'references': ', '.join(set(e['reference'] for e in ph_entries if e['reference']))
                    })
    
    return matched_data

def main():
    """Main function to extract pH data from BRENDA tarball"""
    
    print("BRENDA pH Data Extraction from Tarball")
    print("=" * 50)
    print("Please cite: Schomburg et al., BRENDA in 2021, Nucleic Acids Res. 49(D1)")
    print("=" * 50)
    
    # Look for BRENDA tarball
    tarball_files = list(Path(".").glob("brenda_*.tar.gz"))
    
    if not tarball_files:
        print("\nNo BRENDA tarball found!")
        print("Please download from https://www.brenda-enzymes.org/download.php")
        print("Expected filename pattern: brenda_*.tar.gz")
        return
    
    # Use the first tarball found
    tarball = tarball_files[0]
    print(f"\nFound tarball: {tarball}")
    
    # Extract tarball
    txt_file = extract_brenda_tarball(str(tarball))
    
    if not txt_file:
        print("Failed to extract BRENDA data!")
        return
    
    # Parse BRENDA data
    brenda_data = parse_brenda_flatfile(txt_file)
    
    if not brenda_data:
        print("Failed to parse BRENDA data!")
        return
    
    # Load your PDB structures
    if not os.path.exists("enzyme_ph_dataset.csv"):
        print("\nenzyme_ph_dataset.csv not found!")
        print("Please run the enzyme structure collection script first.")
        return
    
    print("\nLoading PDB structures...")
    pdb_df = pd.read_csv("enzyme_ph_dataset.csv")
    pdb_structures = {
        row['pdb_id']: {
            'ec_number': row['ec_number'],
            'enzyme_name': row.get('enzyme_name', '')
        }
        for _, row in pdb_df.iterrows()
    }
    
    # Match BRENDA data with PDB structures
    print("Matching BRENDA data with PDB structures...")
    matched_data = match_brenda_to_pdb(brenda_data, pdb_structures)
    
    if matched_data:
        df = pd.DataFrame(matched_data)
        output_file = "brenda_pdb_matched.csv"
        df.to_csv(output_file, index=False)
        
        # Create summary statistics
        summary = {
            'total_matches': len(df),
            'unique_ec_numbers': df['ec_number'].nunique(),
            'ph_range': f"{df['ph_optimum'].min():.1f} - {df['ph_optimum'].max():.1f}",
            'avg_ph': df['ph_optimum'].mean(),
            'sources': df['ph_source'].value_counts().to_dict()
        }
        
        print(f"\n{'='*50}")
        print(f"Successfully matched {len(df)} structures with BRENDA pH data")
        print(f"pH range: {summary['ph_range']}")
        print(f"Average pH: {summary['avg_ph']:.1f}")
        print(f"Unique EC numbers: {summary['unique_ec_numbers']}")
        print("\nData sources:")
        for source, count in summary['sources'].items():
            print(f"  {source}: {count}")
        print(f"\nSaved to: {output_file}")
        print(f"{'='*50}")
        
        # Save summary
        with open("brenda_matching_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
    else:
        print("No matches found between BRENDA and PDB structures!")
    
    # Clean up extracted file
    if os.path.exists(txt_file):
        os.remove(txt_file)
        print(f"\nCleaned up temporary file: {txt_file}")

if __name__ == "__main__":
    main()