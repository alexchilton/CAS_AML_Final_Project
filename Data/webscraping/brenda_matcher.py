#!/usr/bin/env python3
"""
Fixed BRENDA matcher that works with your specific CSV format
"""

import tarfile
import pandas as pd
import os
import re
import json

def extract_brenda_tarball(tarball_path="brenda_2024_1.txt.tar.gz"):
    """Extract BRENDA tarball"""
    print(f"Extracting {tarball_path}...")
    
    try:
        with tarfile.open(tarball_path, "r:gz") as tar:
            tar.extractall()
            members = tar.getnames()
            txt_files = [m for m in members if m.endswith('.txt')]
            if txt_files:
                return txt_files[0]
    except Exception as e:
        print(f"Error extracting tarball: {e}")
    
    return None

def parse_brenda_flatfile(filename):
    """Parse BRENDA flatfile format for pH data"""
    print(f"Parsing {filename}...")
    enzymes = {}
    current_ec = None
    
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.rstrip('\n')
            
            # EC number line - handle different possible formats
            if line.startswith("ID\t"):
                parts = line.split('\t')
                if len(parts) > 1:
                    current_ec = parts[1].strip()
                    enzymes[current_ec] = {
                        'name': '',
                        'ph_optimum': [],
                        'ph_range': [],
                        'ph_stability': []
                    }
            
            # Protein name
            elif line.startswith("PR\t") and current_ec:
                parts = line.split('\t', 1)
                if len(parts) > 1:
                    enzymes[current_ec]['name'] = parts[1]
            
            # pH optimum - check multiple possible formats
            elif current_ec and ('PH_OPTIMUM' in line or 'PHO\t' in line):
                ph_data = extract_ph_value(line)
                if ph_data:
                    enzymes[current_ec]['ph_optimum'].extend(ph_data)
            
            # pH range
            elif current_ec and ('PH_RANGE' in line or 'PHR\t' in line):
                ph_data = extract_ph_value(line)
                if ph_data:
                    enzymes[current_ec]['ph_range'].extend(ph_data)
            
            # pH stability
            elif current_ec and ('PH_STABILITY' in line or 'PHS\t' in line):
                ph_data = extract_ph_value(line)
                if ph_data:
                    enzymes[current_ec]['ph_stability'].extend(ph_data)
    
    print(f"Parsed {len(enzymes)} enzymes from BRENDA")
    return enzymes

def extract_ph_value(line):
    """Extract pH values from a line"""
    try:
        # Remove the field identifier
        content = line
        if '\t' in line:
            parts = line.split('\t', 1)
            if len(parts) > 1:
                content = parts[1]
        
        entries = []
        
        # Look for pH values using various patterns
        patterns = [
            r'pH\s*[:=]?\s*(\d+(?:\.\d+)?)',
            r'(\d+(?:\.\d+)?)\s*\(',
            r'\b(\d+(?:\.\d+)?)\b',
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, content)
            for match in matches:
                try:
                    value = float(match)
                    if 0 <= value <= 14:  # Valid pH range
                        entries.append({'value': value})
                except:
                    continue
        
        # Handle pH ranges (e.g., "5.0-7.0")
        range_pattern = r'(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)'
        range_matches = re.findall(range_pattern, content)
        for low, high in range_matches:
            try:
                low_val = float(low)
                high_val = float(high)
                if 0 <= low_val <= 14 and 0 <= high_val <= 14:
                    # Use midpoint of range
                    entries.append({'value': (low_val + high_val) / 2})
            except:
                continue
        
        return entries
    
    except Exception as e:
        return []

def normalize_ec_number(ec):
    """Normalize EC number for matching"""
    # Remove trailing dashes
    ec = ec.rstrip('.-')
    # Remove trailing zeros
    parts = ec.split('.')
    normalized_parts = []
    for part in parts:
        if part.isdigit():
            normalized_parts.append(str(int(part)))
        else:
            normalized_parts.append(part)
    return '.'.join(normalized_parts)

def match_brenda_to_pdb(brenda_data, pdb_df):
    """Match BRENDA pH data with PDB structures"""
    matched_data = []
    unmatched = []
    
    for idx, row in pdb_df.iterrows():
        pdb_id = row['pdb_id']
        ec_number = row['ec_number']
        enzyme_name = row['enzyme_name']
        category = row['category']
        original_ph = row['ph_optimum']
        
        # Try different EC number formats
        ec_variants = [
            ec_number,
            normalize_ec_number(ec_number),
            ec_number.replace('.-', ''),
            ec_number.replace('-', ''),
        ]
        
        matched = False
        for ec_variant in ec_variants:
            if ec_variant in brenda_data:
                brenda_info = brenda_data[ec_variant]
                
                # Get pH values
                ph_values = []
                source_type = None
                
                # Priority: pH optimum > pH range > pH stability
                if brenda_info['ph_optimum']:
                    ph_values = [entry['value'] for entry in brenda_info['ph_optimum'] if 'value' in entry]
                    source_type = 'pH_optimum'
                elif brenda_info['ph_range']:
                    ph_values = [entry['value'] for entry in brenda_info['ph_range'] if 'value' in entry]
                    source_type = 'pH_range'
                elif brenda_info['ph_stability']:
                    ph_values = [entry['value'] for entry in brenda_info['ph_stability'] if 'value' in entry]
                    source_type = 'pH_stability'
                
                if ph_values:
                    # Use median pH value
                    median_ph = sorted(ph_values)[len(ph_values)//2]
                    
                    matched_data.append({
                        'category': category,
                        'enzyme_name': enzyme_name,
                        'pdb_id': pdb_id,
                        'ec_number': ec_number,
                        'original_ph': original_ph,
                        'brenda_ph': median_ph,
                        'ph_source': f'BRENDA_{source_type}',
                        'ph_values_count': len(ph_values),
                        'brenda_ec': ec_variant
                    })
                    matched = True
                    break
        
        if not matched:
            unmatched.append({
                'enzyme_name': enzyme_name,
                'ec_number': ec_number,
                'pdb_id': pdb_id
            })
    
    return matched_data, unmatched

def main():
    """Main function to extract pH data from BRENDA"""
    
    print("BRENDA pH Data Extraction")
    print("=" * 50)
    
    # Check for required files
    #filename='expanded_enzyme_dataset.csv'
    if not os.path.exists(filename):     ## Set here the name of the csv file
        print("ERROR: enzyme_ph_dataset.csv not found!")
        return
    
    # Find BRENDA tarball
    tarball_files = [f for f in os.listdir('.') if f.endswith('.tar.gz') and 'brenda' in f.lower()]
    if not tarball_files:
        print("ERROR: No BRENDA tarball found!")
        return
    
    tarball = tarball_files[0]
    print(f"Found tarball: {tarball}")
    
    # Extract tarball
    txt_file = extract_brenda_tarball(tarball)
    if not txt_file:
        print("Failed to extract BRENDA data!")
        return
    
    # Parse BRENDA data
    brenda_data = parse_brenda_flatfile(txt_file)
    
    # Load PDB structures
    print("\nLoading PDB structures...")
    pdb_df = pd.read_csv(filename)
    print(f"Loaded {len(pdb_df)} structures")
    
    # Match BRENDA data with PDB structures
    print("\nMatching BRENDA data with PDB structures...")
    matched_data, unmatched_data = match_brenda_to_pdb(brenda_data, pdb_df)
    
    if matched_data:
        # Create results DataFrame
        results_df = pd.DataFrame(matched_data)
        
        # Save matched data
        output_file = "brenda_pdb_matched.csv"
        results_df.to_csv(output_file, index=False)
        
        # Create enhanced dataset
        enhanced_df = pdb_df.copy()
        for idx, row in results_df.iterrows():
            mask = enhanced_df['pdb_id'] == row['pdb_id']
            enhanced_df.loc[mask, 'ph_optimum'] = row['brenda_ph']
            enhanced_df.loc[mask, 'ph_source'] = row['ph_source']
        
        enhanced_df.to_csv(output_name, index=False)
        
        # Create summary
        summary = {
            'total_matches': len(results_df),
            'unmatched': len(unmatched_data),
            'match_rate': f"{len(results_df)/len(pdb_df)*100:.1f}%",
            'ph_range': f"{results_df['brenda_ph'].min():.1f} - {results_df['brenda_ph'].max():.1f}",
            'avg_ph_change': f"{abs(results_df['original_ph'] - results_df['brenda_ph']).mean():.2f}",
            'categories_matched': results_df['category'].nunique(),
            'sources': results_df['ph_source'].value_counts().to_dict()
        }
        
        print(f"\n{'='*50}")
        print(f"Successfully matched {len(results_df)} structures with BRENDA")
        print(f"Match rate: {summary['match_rate']}")
        print(f"pH range from BRENDA: {summary['ph_range']}")
        print(f"Average pH change: {summary['avg_ph_change']}")
        print(f"\nFiles created:")
        print(f"  - {output_file} (detailed matches)")
        print(f"  - enzyme_ph_dataset_enhanced.csv (updated dataset)")
        print(f"{'='*50}")
        
        # Save summary
        with open("brenda_matching_summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        
        # Report unmatched
        if unmatched_data:
            unmatched_df = pd.DataFrame(unmatched_data)
            unmatched_df.to_csv("unmatched_enzymes.csv", index=False)
            print(f"\nUnmatched enzymes saved to: unmatched_enzymes.csv")
            print(f"Example unmatched: {unmatched_data[:3]}")
    else:
        print("\nNo matches found!")
        print("Debugging: First few EC numbers from your dataset:")
        for ec in pdb_df['ec_number'].unique()[:5]:
            print(f"  {ec}")
        print("\nFirst few EC numbers from BRENDA:")
        for ec in list(brenda_data.keys())[:5]:
            print(f"  {ec}")
    
    # Clean up
    if os.path.exists(txt_file):
        os.remove(txt_file)
        print(f"\nCleaned up temporary file: {txt_file}")

if __name__ == "__main__":
    filename='expanded_enzyme_dataset.csv'
    output_name= 'matched_pdb_brenda_files.csv'
    main()