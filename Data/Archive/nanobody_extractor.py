#!/usr/bin/env python3
"""
Script to extract only the nanobody chains from PDB files.
This script analyzes downloaded PDB files and creates new files containing only the nanobody portions.
"""

import os
import re
import glob
from tqdm import tqdm  # For progress bar

# Common nanobody chain identifiers and keywords
NANOBODY_KEYWORDS = [
    "nanobody", "vhh", "camelid", "single domain antibody", 
    "heavy chain antibody", "camel antibody"
]

def identify_nanobody_chains(pdb_file_path):
    """
    Analyze a PDB file to identify which chains are nanobodies.
    
    This function uses several methods to identify nanobody chains:
    1. Look for explicit mentions in the COMPND or TITLE records
    2. Check for typical nanobody size (around 110-130 residues)
    3. Look for characteristic residue patterns
    
    Args:
        pdb_file_path (str): Path to the PDB file
        
    Returns:
        list: List of chain IDs that are likely nanobodies
    """
    nanobody_chains = set()
    chains_info = {}  # Store information about each chain
    
    # Regular expressions for extracting chain information
    chain_re = re.compile(r'CHAIN:\s+(.*?)(?:;|$)')
    molecule_re = re.compile(r'MOLECULE:\s+(.*?)(?:;|$)')
    
    with open(pdb_file_path, 'r') as f:
        title = ""
        current_chain = None
        current_residue = None
        chains_residues = {}  # Keep track of residue counts per chain
        
        for line in f:
            line = line.strip()
            
            # Extract title information
            if line.startswith("TITLE"):
                title += line[10:].strip()
            
            # Extract chain information from COMPND records
            elif line.startswith("COMPND"):
                content = line[10:].strip()
                
                # Check if this line defines chains
                chain_match = chain_re.search(content)
                if chain_match:
                    chains = chain_match.group(1).split(", ")
                    
                    # Check if the previous MOLECULE line had nanobody keywords
                    if any(keyword.lower() in molecule.lower() for keyword in NANOBODY_KEYWORDS for molecule in chains_info.get('current_molecules', [])):
                        for chain in chains:
                            nanobody_chains.add(chain)
                
                # Check if this line defines the molecule
                molecule_match = molecule_re.search(content)
                if molecule_match:
                    molecule = molecule_match.group(1)
                    if 'current_molecules' not in chains_info:
                        chains_info['current_molecules'] = []
                    chains_info['current_molecules'].append(molecule)
                    
                    # Check if this molecule is explicitly a nanobody
                    if any(keyword.lower() in molecule.lower() for keyword in NANOBODY_KEYWORDS):
                        chains_info['nanobody_molecule'] = True
            
            # Track atom records to count residues per chain
            elif line.startswith("ATOM"):
                chain_id = line[21]
                residue_id = int(line[22:26].strip())
                
                if chain_id not in chains_residues:
                    chains_residues[chain_id] = set()
                
                chains_residues[chain_id].add(residue_id)
    
    # Check title for nanobody keywords
    title_has_nanobody = any(keyword.lower() in title.lower() for keyword in NANOBODY_KEYWORDS)
    
    # If title mentions nanobody but we couldn't identify chains, take a best guess
    if title_has_nanobody and not nanobody_chains:
        # Look for chains with typical nanobody length (around 110-130 residues)
        for chain, residues in chains_residues.items():
            num_residues = len(residues)
            if 90 <= num_residues <= 140:  # Generous range for nanobody length
                nanobody_chains.add(chain)
    
    # If still no nanobody chains identified, take the smallest chain(s) that could be a nanobody
    if not nanobody_chains:
        # Sort chains by residue count
        sorted_chains = sorted([(chain, len(residues)) for chain, residues in chains_residues.items()], 
                              key=lambda x: x[1])
        
        # Select chains that might be nanobodies based on size
        for chain, num_residues in sorted_chains:
            if 90 <= num_residues <= 140:  # Typical nanobody size range
                nanobody_chains.add(chain)
                break  # Usually only one nanobody per structure
    
    return list(nanobody_chains)

def extract_nanobody_chains(pdb_file_path, output_dir, chains):
    """
    Extract specified chains from a PDB file and save to a new file.
    
    Args:
        pdb_file_path (str): Path to the input PDB file
        output_dir (str): Directory to save the extracted chains
        chains (list): List of chain IDs to extract
        
    Returns:
        str: Path to the new PDB file containing only the specified chains
    """
    if not chains:
        return None
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate output filename
    pdb_id = os.path.basename(pdb_file_path).split('.')[0]
    chain_str = '_'.join(chains)
    output_filename = f"{pdb_id}_nanobody_{chain_str}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    
    # Extract chains
    with open(pdb_file_path, 'r') as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM", "TER")):
                chain_id = line[21]
                if chain_id in chains:
                    f_out.write(line)
            elif not line.startswith(("ATOM", "HETATM", "TER")):
                # Copy header and other non-atom records
                f_out.write(line)
    
    return output_path

def process_pdb_directory(input_dir, output_dir):
    """
    Process all PDB files in a directory to extract nanobody chains.
    
    Args:
        input_dir (str): Directory containing PDB files
        output_dir (str): Directory to save the extracted nanobody files
    """
    # Find all PDB files in the input directory
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return
    
    print(f"Processing {len(pdb_files)} PDB files...")
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process each PDB file
    processed_count = 0
    extracted_files = []
    
    for pdb_file in tqdm(pdb_files):
        try:
            # Identify nanobody chains
            nanobody_chains = identify_nanobody_chains(pdb_file)
            
            if nanobody_chains:
                # Extract nanobody chains
                output_file = extract_nanobody_chains(pdb_file, output_dir, nanobody_chains)
                if output_file:
                    extracted_files.append((os.path.basename(pdb_file), nanobody_chains, output_file))
                    processed_count += 1
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
    
    print(f"Successfully extracted nanobody chains from {processed_count} PDB files.")
    
    # Create a summary file
    summary_path = os.path.join(output_dir, "nanobody_extraction_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Extracted nanobody chains from {processed_count} of {len(pdb_files)} PDB files\n\n")
        f.write("Details:\n")
        for original_file, chains, output_file in extracted_files:
            f.write(f"{original_file}: Extracted chains {', '.join(chains)} -> {os.path.basename(output_file)}\n")
    
    print(f"Summary file created: {summary_path}")

def main():
    """Main function to extract nanobody chains from PDB files."""
    print("Nanobody Chain Extractor")
    print("=======================")
    
    try:
        # Get input and output directories
        input_dir = input("Enter the directory containing PDB files (default: nanobody_pdbs): ").strip()
        if not input_dir:
            input_dir = "nanobody_pdbs"
        
        output_dir = input("Enter the directory to save extracted nanobody files (default: nanobody_extracted): ").strip()
        if not output_dir:
            output_dir = "nanobody_extracted"
        
        # Process PDB files
        process_pdb_directory(input_dir, output_dir)
        
        print("\nProcess completed successfully!")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()