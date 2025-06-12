from Bio import PDB
from Bio.PDB.DSSP import DSSP
import networkx as nx

def calculate_secondary_structure(pdb_path):
    """
    Calculate secondary structure using DSSP.
    
    Parameters:
    pdb_path -- Path to the PDB file
    
    Returns:
    Dictionary mapping (chain_id, residue_id) to secondary structure
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        model = structure[0]
        
        try:
            # Try direct DSSP call
            dssp = DSSP(model, pdb_path)
        except Exception as e1:
            try:
                # Fallback: Try with DSSP executable if installed
                from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                dssp_dict, dssp_keys = dssp_dict_from_pdb_file(pdb_path)
                
                # Convert to format similar to DSSP object
                dssp = {}
                for key in dssp_keys:
                    dssp[key] = dssp_dict[key]
            except Exception as e2:
                print(f"DSSP calculation failed: {str(e1)} | {str(e2)}")
                # If all methods fail, use a simple fallback
                return {}
        
        ss_data = {}
        
        for key in dssp.keys():
            if isinstance(key, tuple):
                # Newer DSSP version
                chain_id = key[0]
                res_id = key[1][1]  # Residue number
                if isinstance(dssp[key], dict):
                    ss = dssp[key].get('secstruct', '?')
                else:
                    # Older DSSP version
                    ss = dssp[key][2]  # Secondary structure
            else:
                # Alternative key format
                chain_id, res_id = key[0], key[1]
                ss = dssp[key][2]
                
            ss_data[(chain_id, res_id)] = ss
        
        return ss_data
    except Exception as e:
        print(f"Secondary structure calculation failed: {str(e)}")
        return {}

def assign_charges(chain_atoms):
    """Assign partial charges to atoms."""
    charges = {}
    
    # Simplified charge assignment based on atom types
    charge_map = {
        'N': -0.5,   # Backbone nitrogen
        'C': 0.5,    # Backbone carbon
        'O': -0.5,   # Carbonyl oxygen
        'S': 0.0,    # Sulfur
        'P': 0.5     # Phosphorus
    }
    
    for i, atom in enumerate(chain_atoms):
        chain_id, res_id, res_name, atom_name, element, position = atom
        charges[i] = charge_map.get(element, 0.0)
    
    return charges

def identify_hydrophobic_residues(chain_residues):
    """Identify hydrophobic residues in a chain."""
    hydrophobic_aas = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
    hydrophobic_res_ids = set()
    
    for res_id, res_name in chain_residues:
        if res_name in hydrophobic_aas:
            hydrophobic_res_ids.add(res_id)
    
    return hydrophobic_res_ids