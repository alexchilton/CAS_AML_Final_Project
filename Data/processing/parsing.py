import gc
from Bio import PDB

def parse_pdb(pdb_filename):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("Protein", pdb_filename)
    
    atoms = []
    residues = []
    
    # Keep track of processed residues by chain
    all_residues_by_chain = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            chain_residues = []
            
            for residue in chain:
                # Get residue details
                residue_id = residue.get_id()[1]  # The sequence number
                residue_name = residue.get_resname()
                
                # Add residue to this chain's list
                chain_residues.append((residue_id, residue_name))
                
                # Process all atoms in this residue
                for atom in residue:
                    atom_name = atom.get_name()
                    element = atom.element
                    position = atom.get_coord()
                    # Include chain ID in atom information for clarity
                    atoms.append((chain_id, residue_id, residue_name, atom_name, element, position))
            
            # Store all residues for this chain
            all_residues_by_chain[chain_id] = chain_residues
    
    # Option 1: Return separate residue lists for each chain --> better for our pipeline
    return atoms, all_residues_by_chain

def parse_basic_structure(pdb_path):
    """
    Parse basic PDB structure with low memory usage.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
        
    Returns:
    --------
    tuple
        (structure object, residues by chain)
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("temp", pdb_path)
        
        # Extract residues by chain
        residues_by_chain = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                residues_by_chain[chain_id] = []
                
                for residue in chain:
                    # Skip non-standard residues
                    if residue.get_id()[0] != " ":
                        continue
                    
                    res_id = residue.get_id()[1]
                    res_name = residue.get_resname()
                    residues_by_chain[chain_id].append((res_id, res_name))
        
        return structure, residues_by_chain
    
    except Exception as e:
        print(f"Error parsing structure: {str(e)}")
        return None, None


def extract_chain_atoms(pdb_path, chain_id):
    """Extract atoms for a specific chain only."""
    atoms = []
    
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("temp", pdb_path)
    
    for model in structure:
        for chain in model:
            if chain.get_id() != chain_id:
                continue
                
            for residue in chain:
                # Skip non-standard residues
                if residue.get_id()[0] != " ":
                    continue
                    
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()
                
                for atom in residue:
                    atom_name = atom.get_name()
                    element = atom.element
                    position = atom.get_coord()
                    
                    atoms.append((chain_id, res_id, res_name, atom_name, element, position))
    
    # Clear reference to structure
    del structure
    gc.collect()
    
    return atoms

