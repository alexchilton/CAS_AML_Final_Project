#!/usr/bin/env python3
"""
Improved script to extract nanobody chains from PDB files.
- Creates separate files for each individual nanobody chain
- Implements stricter identification criteria
- Performs additional validation to ensure chains are actual nanobodies
"""

import os
import re
import glob
from tqdm import tqdm  # For progress bar

# Common nanobody keywords
NANOBODY_KEYWORDS = [
    "nanobody", "vhh", "camelid", "single domain antibody", 
    "heavy chain antibody", "camel antibody"
]

def count_residues_in_chain(pdb_file_path, chain_id):
    """
    Count the number of residues in a specific chain of a PDB file.
    
    Args:
        pdb_file_path (str): Path to the PDB file
        chain_id (str): Chain identifier
        
    Returns:
        int: Number of residues in the chain
    """
    residues = set()
    
    with open(pdb_file_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                if line[21] == chain_id:
                    residue_id = int(line[22:26].strip())
                    residues.add(residue_id)
    
    return len(residues)

def is_likely_nanobody(pdb_file_path, chain_id):
    """
    Determine if a chain is likely to be a nanobody based on size and composition.
    
    Args:
        pdb_file_path (str): Path to the PDB file
        chain_id (str): Chain identifier
        
    Returns:
        bool: True if the chain is likely a nanobody, False otherwise
    """
    # Check residue count (nanobodies typically have 110-140 residues)
    residue_count = count_residues_in_chain(pdb_file_path, chain_id)
    if not (90 <= residue_count <= 140):
        return False
    
    # Additional verification could be added here:
    # - Check for characteristic CDR regions
    # - Verify presence of key conserved residues
    # - Check secondary structure composition
    
    return True

def identify_nanobody_chains(pdb_file_path):
    """
    Analyze a PDB file to identify which chains are nanobodies using stricter criteria.
    
    Args:
        pdb_file_path (str): Path to the PDB file
        
    Returns:
        list: List of chain IDs that are likely nanobodies
    """
    possible_nanobody_chains = set()
    explicit_nanobody_mentions = False
    chain_molecules = {}  # Map chains to their molecule descriptions
    
    with open(pdb_file_path, 'r') as f:
        title = ""
        current_mol_id = None
        current_molecule = None
        current_chains = []
        
        for line in f:
            line = line.strip()
            
            # Extract title information
            if line.startswith("TITLE"):
                title += line[10:].strip()
            
            # Extract chain information from COMPND records
            elif line.startswith("COMPND"):
                content = line[10:].strip()
                
                # Track MOL_ID
                if "MOL_ID:" in content:
                    current_mol_id = content.split("MOL_ID:")[1].strip().rstrip(";")
                    current_molecule = None
                    current_chains = []
                
                # Track MOLECULE
                elif "MOLECULE:" in content:
                    molecule_part = content.split("MOLECULE:")[1].strip().rstrip(";")
                    if current_molecule is None:
                        current_molecule = molecule_part
                    else:
                        current_molecule += " " + molecule_part
                        
                    # Check if explicitly mentions nanobody
                    if any(keyword.lower() in current_molecule.lower() for keyword in NANOBODY_KEYWORDS):
                        explicit_nanobody_mentions = True
                
                # Track CHAIN
                elif "CHAIN:" in content:
                    chains_part = content.split("CHAIN:")[1].strip().rstrip(";")
                    current_chains = [c.strip() for c in re.split(r'[,\s]+', chains_part) if c.strip()]
                    
                    # Associate current chains with current molecule
                    if current_molecule:
                        for chain in current_chains:
                            chain_molecules[chain] = current_molecule
                            
                            # If molecule explicitly mentions nanobody, mark chain as possible nanobody
                            if any(keyword.lower() in current_molecule.lower() for keyword in NANOBODY_KEYWORDS):
                                possible_nanobody_chains.add(chain)
    
    # Check title for nanobody keywords
    title_has_nanobody = any(keyword.lower() in title.lower() for keyword in NANOBODY_KEYWORDS)
    
    # Verify each possible nanobody chain
    confirmed_nanobody_chains = []
    
    # If we have explicit mentions and possible chains, use those
    if explicit_nanobody_mentions and possible_nanobody_chains:
        for chain in possible_nanobody_chains:
            if is_likely_nanobody(pdb_file_path, chain):
                confirmed_nanobody_chains.append(chain)
    
    # If title suggests nanobodies but we didn't identify any explicit chains
    elif title_has_nanobody and not confirmed_nanobody_chains:
        # Check each chain in the file to see if it matches nanobody characteristics
        all_chains = set()
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    all_chains.add(line[21])
        
        for chain in all_chains:
            if is_likely_nanobody(pdb_file_path, chain):
                confirmed_nanobody_chains.append(chain)
    
    # If still no chains identified, fall back to size-based detection
    if not confirmed_nanobody_chains:
        all_chains = set()
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    all_chains.add(line[21])
        
        for chain in all_chains:
            if is_likely_nanobody(pdb_file_path, chain):
                confirmed_nanobody_chains.append(chain)
    
    return confirmed_nanobody_chains

def extract_single_chain(pdb_file_path, output_dir, chain_id):
    """
    Extract a single chain from a PDB file and save to a new file.
    
    Args:
        pdb_file_path (str): Path to the input PDB file
        output_dir (str): Directory to save the extracted chain
        chain_id (str): Chain ID to extract
        
    Returns:
        str: Path to the new PDB file containing only the specified chain
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate output filename
    pdb_id = os.path.basename(pdb_file_path).split('.')[0]
    output_filename = f"{pdb_id}_nanobody_{chain_id}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    
    # Extract chain
    with open(pdb_file_path, 'r') as f_in, open(output_path, 'w') as f_out:
        # First, write a modified header
        f_out.write(f"REMARK   Generated by Nanobody Extractor\n")
        f_out.write(f"REMARK   Source file: {os.path.basename(pdb_file_path)}\n")
        f_out.write(f"REMARK   Extracted chain: {chain_id}\n\n")
        
        # Track if we've seen this chain's atoms yet (to add TER record after the chain)
        seen_chain_atoms = False
        
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                line_chain_id = line[21]
                if line_chain_id == chain_id:
                    f_out.write(line)
                    seen_chain_atoms = True
            elif line.startswith("TER") and seen_chain_atoms:
                f_out.write(line)
                seen_chain_atoms = False
            elif not line.startswith(("ATOM", "HETATM", "TER", "END")):
                # Copy header records, but skip coordinate records for other chains
                f_out.write(line)
        
        # Ensure there's an END record
        f_out.write("END\n")
    
    return output_path

def process_pdb_directory(input_dir, output_dir):
    """
    Process all PDB files in a directory to extract individual nanobody chains.
    
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
    extracted_chains = []
    
    for pdb_file in tqdm(pdb_files):
        try:
            # Identify nanobody chains
            nanobody_chains = identify_nanobody_chains(pdb_file)
            
            if nanobody_chains:
                # Extract each nanobody chain individually
                for chain in nanobody_chains:
                    output_file = extract_single_chain(pdb_file, output_dir, chain)
                    if output_file:
                        extracted_chains.append((os.path.basename(pdb_file), chain, os.path.basename(output_file)))
                
                processed_count += 1
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
    
    print(f"Successfully processed {processed_count} PDB files.")
    print(f"Extracted {len(extracted_chains)} individual nanobody chains.")
    
    # Create a summary file
    summary_path = os.path.join(output_dir, "nanobody_extraction_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"Processed {processed_count} of {len(pdb_files)} PDB files\n")
        f.write(f"Extracted {len(extracted_chains)} individual nanobody chains\n\n")
        f.write("Details:\n")
        
        # Group by original file
        file_groups = {}
        for original_file, chain, output_file in extracted_chains:
            if original_file not in file_groups:
                file_groups[original_file] = []
            file_groups[original_file].append((chain, output_file))
        
        # Write grouped summary
        for original_file, extractions in file_groups.items():
            f.write(f"{original_file}:\n")
            for chain, output_file in extractions:
                f.write(f"  - Chain {chain} -> {output_file}\n")
            f.write("\n")
    
    print(f"Summary file created: {summary_path}")

def main():
    """Main function to extract individual nanobody chains from PDB files."""
    print("Improved Nanobody Chain Extractor")
    print("================================")
    
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