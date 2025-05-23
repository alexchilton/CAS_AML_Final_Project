# Protein Graph Data Format Documentation

## Overview
This document describes the data structures produced by the pdb processing pipeline. The output consists of NetworkX graphs and associated physical properties stored in pickled dictionaries.

## File Structure
For each PDB file processed, the following files are generated:

### Main Files
- `{pdb_id}_graph.pkl`: NetworkX graph representing residue-level connectivity for the entire protein
- `{pdb_id}_data.pkl`: Main dictionary containing aggregated data for all chains

### Chain-Specific Files
For each chain (e.g., H, L) in the protein, separate files are generated to optimize memory usage:

- `{pdb_id}_{chain_id}_atoms.pkl`: List of all atoms in this chain
- `{pdb_id}_{chain_id}_backbone.pkl`: Dictionary of backbone atom coordinates
- `{pdb_id}_{chain_id}_ss.pkl`: Secondary structure assignments
- `{pdb_id}_{chain_id}_bond_lengths.pkl`: Dictionary of bond lengths
- `{pdb_id}_{chain_id}_angles.pkl`: Dictionary of bond angles
- `{pdb_id}_{chain_id}_torsions.pkl`: Dictionary of torsion angles (backbone only)
- `{pdb_id}_{chain_id}_charges.pkl`: Dictionary of atomic partial charges
- `{pdb_id}_{chain_id}_hydrophobic.pkl`: Set of hydrophobic residue IDs
- `{pdb_id}_{chain_id}_edge_pairs.pkl`: Atom connectivity for this chain
- `{pdb_id}_{chain_id}_index.json`: Index file linking all chain-specific files


## NetworkX Graph Structure

### Nodes
Each node in the graph represents a single amino acid residue with the following attributes:

- **Node ID format**: `"{chain_id}:{residue_name}:{residue_number}"` (e.g., "A:VAL:23")
- **Attributes**:
  - `chain_id`: Chain identifier (e.g., "A", "B")
  - `residue_number`: Position in the protein sequence
  - `residue_name`: Three-letter amino acid code
  - `ss`: Secondary structure assignment (DSSP code)
    - H: Alpha helix
    - E: Extended strand (beta sheet)
    - B: Beta bridge
    - G: 3-10 helix
    - T: Turn
    - S: Bend
    - ?: Unknown
  - `has_backbone`: Boolean indicating if backbone atoms are present
  - `coords`: 3D coordinates (typically CA atom position)
  - Various atom coordinates (when available): `CA`, `N`, `C`, `O`

### Edges
Edges represent connections between residues:

- **Types**:
  - Sequential connections (peptide bonds) between consecutive residues
  - Spatial proximity connections (residues with CA atoms within 7Å)
- **Attributes**:
  - `kind`: Edge type (e.g., `{"peptide_bond"}` or `{"contact"}`)
  - `distance`: For spatial proximity edges

## Data Dictionary Structure

The `{pdb_id}_data.pkl` file contains a dictionary with the following structure:

```python
{
    'pdb_id': str,  # PDB identifier
    
    'atoms': list,  # List of atom tuples (chain_id, res_id, res_name, atom_name, element, [x,y,z])
    
    'residues_by_chain': {
        'chain_id': [(res_id, res_name), ...],
        # Dictionary mapping chain IDs to lists of residue information
    },
    
    'secondary_structure': {
        'chain_id': {res_id: ss_code, ...},
        # Dictionary of secondary structure assignments
    },
    
    'backbone_atoms': {
        'chain_id': {
            res_id: {
                'N': [x, y, z],
                'CA': [x, y, z],
                'C': [x, y, z],
                'O': [x, y, z]
            },
            # Dictionary mapping residue IDs to backbone atom coordinates
        }
    },
    
    'edge_indices': {
        'chain_id': array(shape=(2, num_edges)),
        # Atom connectivity for each chain
    },
    
    'bond_lengths': {
        'chain_id': {
            (atom_i, atom_j): length,
            # Dictionary mapping atom pairs to bond lengths
        }
    },
    
    'angles': {
        'chain_id': {
            (atom_i, atom_j, atom_k): angle_in_degrees,
            # Dictionary mapping atom triplets to bond angles
        }
    },
    
    'torsions': {
        'chain_id': {
            (atom_i, atom_j, atom_k, atom_l): torsion_in_degrees,
            # Dictionary mapping atom quadruplets to torsion angles
        }
    },
    
    'charges': {
        'chain_id': {
            atom_index: charge_value,
            # Dictionary mapping atom indices to partial charges
        }
    },
    
    'hydrophobic_residues': {
        'chain_id': {res_id, res_id, ...},
        # Set of hydrophobic residue IDs per chain
    }
}

```

## Physical Properties

- **Bond lengths**: Euclidean distances between connected atoms
- **Angles**: Angles in degrees between triplets of connected atoms
- **Torsions**: Dihedral angles in degrees (primarily backbone phi/psi/omega angles)
- **Charges**: Simplified partial charges assigned based on atom type
- **Hydrophobic residues**: Residues identified as hydrophobic (ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR)

## Memory Optimization

The processing pipeline implements memory-efficient strategies:

1. Chain-by-chain processing
2. Sparse representation of physical properties
3. Distance cutoffs for interactions (typically 2-5Å)
4. Selective calculation of torsion angles (backbone only)


