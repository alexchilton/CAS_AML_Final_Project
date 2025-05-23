import os
import numpy as np
import pickle
import json
import traceback
import gc
from tqdm import tqdm
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import networkx as nx

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

def process_pdb_full_features_memory_efficient(folder_path, output_path=None, max_file_size_mb=25):
    """
    Process PDB files with all features while being memory-efficient.
    
    Parameters:
    -----------
    folder_path : str
        Path to folder with PDB files
    output_path : str, optional
        Output directory
    max_file_size_mb : float
        Skip files larger than this
        
    Returns:
    --------
    dict
        Processing summary
    """
    if output_path is None:
        output_path = os.path.join(folder_path, "processed_data")
    
    os.makedirs(output_path, exist_ok=True)
    
    # Get PDB files
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]
    print(f"Found {len(pdb_files)} PDB files")
    
    summary = {}
    
    for pdb_file in pdb_files:
        pdb_path = os.path.join(folder_path, pdb_file)
        pdb_id = os.path.splitext(pdb_file)[0]
        
        # Create PDB-specific output directory
        pdb_output = os.path.join(output_path, pdb_id)
        os.makedirs(pdb_output, exist_ok=True)
        
        try:
            # Check file size
            file_size_mb = os.path.getsize(pdb_path) / (1024 * 1024)
            if file_size_mb > max_file_size_mb:
                print(f"Skipping {pdb_file} - too large ({file_size_mb:.2f} MB)")
                summary[pdb_id] = {"status": "skipped", "size_mb": file_size_mb}
                continue
            
            print(f"Processing {pdb_file} ({file_size_mb:.2f} MB)")
            
            # Step 1: Parse basic structure and create graph (low memory impact)
            structure, residues_by_chain = parse_basic_structure(pdb_path)
            
            if not structure or not residues_by_chain:
                print(f"Failed to parse {pdb_file}")
                summary[pdb_id] = {"status": "error", "error": "Failed to parse structure"}
                continue
            
            print(f"  Found {len(residues_by_chain)} chains")
            
            # Step 2: Calculate secondary structure once (low memory impact)
            ss_data = calculate_secondary_structure(pdb_path)
            print(f"  Calculated secondary structure for {len(ss_data)} residues")
            
            # Step 3: Create protein graph (moderate memory impact)
            nx_graph = create_protein_graph(pdb_path)
            nx_graph = add_structure_to_graph(nx_graph, ss_data)
            print(f"  Created graph with {nx_graph.number_of_nodes()} nodes and {nx_graph.number_of_edges()} edges")
            
            # Save graph immediately to free memory
            graph_path = os.path.join(pdb_output, f"{pdb_id}_graph.pkl")
            with open(graph_path, 'wb') as f:
                pickle.dump(nx_graph, f)
            print(f"  Saved graph to {graph_path}")
            
            # Create PDB data dictionary (without the graph)
            pdb_data = {
                'pdb_id': pdb_id,
                'residues_by_chain': residues_by_chain,
                'secondary_structure': {},
                'backbone_atoms': {},
                'edge_indices': {},
                'bond_lengths': {},
                'angles': {},
                'torsions': {},
                'charges': {},
                'hydrophobic_residues': {}
            }
            
            # Process each chain separately with incremental save
            for chain_id in residues_by_chain:
                print(f"  Processing chain {chain_id}...")
                
                # Create output path for this chain
                chain_output = os.path.join(pdb_output, f"{pdb_id}_{chain_id}")
                
                # Step 4: Process one chain at a time with memory-efficient functions
                #         Save results incrementally to disk
                process_single_chain_full_features(
                    pdb_path, 
                    chain_id, 
                    chain_output, 
                    ss_data,
                    incremental_save=True
                )
                
                # After processing is complete, load the chain results
                chain_data = load_chain_data(chain_output)
                
                # Add to main data dictionary
                if chain_data:
                    for key in ['secondary_structure', 'backbone_atoms', 'edge_indices', 
                               'bond_lengths', 'angles', 'torsions', 'charges', 
                               'hydrophobic_residues']:
                        if key in chain_data:
                            pdb_data[key][chain_id] = chain_data[key]
                
                # Force garbage collection
                gc.collect()
            
            # Save the main data dictionary
            data_path = os.path.join(pdb_output, f"{pdb_id}_data.pkl")
            with open(data_path, 'wb') as f:
                pickle.dump(pdb_data, f)
            print(f"  Saved data to {data_path}")
            
            # Update summary
            summary[pdb_id] = {
                "status": "success",
                "chains": list(residues_by_chain.keys()),
                "size_mb": file_size_mb,
                "num_nodes": nx_graph.number_of_nodes(),
                "num_edges": nx_graph.number_of_edges()
            }
            
            # Clear memory for next PDB
            del nx_graph
            del pdb_data
            gc.collect()
            
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
            traceback.print_exc()
            summary[pdb_id] = {"status": "error", "error": str(e)}
    
    # Save summary
    summary_path = os.path.join(output_path, "processing_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    return summary

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

def process_single_chain_full_features(pdb_path, chain_id, output_prefix, ss_data, incremental_save=True):
    """
    Process a single chain with all features but in a memory-efficient way.
    
    Parameters:
    -----------
    pdb_path : str
        Path to PDB file
    chain_id : str
        Chain ID to process
    output_prefix : str
        Prefix for output files
    ss_data : dict
        Secondary structure data
    incremental_save : bool
        Whether to save results incrementally
        
    Returns:
    --------
    dict
        Chain data if not saved incrementally
    """
    try:
        # Step 1: Extract atoms for this chain only
        atoms = extract_chain_atoms(pdb_path, chain_id)
        
        if not atoms:
            print(f"  No atoms found for chain {chain_id}")
            return {} if not incremental_save else None
        
        print(f"  Extracted {len(atoms)} atoms for chain {chain_id}")
        
        # Save atoms if incremental
        if incremental_save:
            with open(f"{output_prefix}_atoms.pkl", 'wb') as f:
                pickle.dump(atoms, f)
        
        # Step 2: Get chain residues
        chain_residues = []
        res_ids_seen = set()
        
        for atom in atoms:
            _, res_id, res_name, _, _, _ = atom
            if res_id not in res_ids_seen:
                chain_residues.append((res_id, res_name))
                res_ids_seen.add(res_id)
        
        # Step 3: Extract backbone atoms (low memory impact)
        backbone = extract_backbone_atoms(atoms)
        print(f"  Found backbone atoms for {len(backbone)} residues")
        
        if incremental_save:
            with open(f"{output_prefix}_backbone.pkl", 'wb') as f:
                pickle.dump(backbone, f)
        
        # Step 4: Get secondary structure for this chain
        chain_ss = {}
        for res_id, _ in chain_residues:
            chain_ss[res_id] = ss_data.get((chain_id, res_id), '?')
        
        if incremental_save:
            with open(f"{output_prefix}_ss.pkl", 'wb') as f:
                pickle.dump(chain_ss, f)
        
        # Step 5: Extract atomic positions (needed for multiple calculations)
        positions = np.array([atom[5] for atom in atoms])
        
        # Step 6: Create edge indices for different calculations
        print("  Creating edge indices...")
        
        # Pairs for bond lengths (with efficient memory)
        edge_index_pairs = create_edge_index_memory_efficient(atoms, positions, mode='pairs')
        
        if incremental_save:
            with open(f"{output_prefix}_edge_pairs.pkl", 'wb') as f:
                pickle.dump(edge_index_pairs, f)
        
        # Step 7: Calculate bond lengths (moderate memory impact)
        print("  Calculating bond lengths...")
        bond_lengths = calculate_bond_lengths_efficient(atoms, edge_index_pairs)
        
        if incremental_save:
            with open(f"{output_prefix}_bond_lengths.pkl", 'wb') as f:
                pickle.dump(bond_lengths, f)
            
            # Clear memory
            del edge_index_pairs
            gc.collect()
        
        # Step 8: Calculate angles (higher memory impact)
        print("  Calculating angles...")
        angles = calculate_angles_memory_efficient(atoms, positions, backbone_only=True)
        
        if incremental_save:
            with open(f"{output_prefix}_angles.pkl", 'wb') as f:
                pickle.dump(angles, f)
            
            # Clear memory
            gc.collect()
        
        # Step 9: Calculate torsions (higher memory impact)
        print("  Calculating torsions...")
        torsions = calculate_torsions_memory_efficient(atoms, positions)
        
        if incremental_save:
            with open(f"{output_prefix}_torsions.pkl", 'wb') as f:
                pickle.dump(torsions, f)
            
            # Clear memory
            gc.collect()
        
        # Step 10: Calculate charges and hydrophobicity (low memory impact)
        print("  Assigning charges and identifying hydrophobic residues...")
        charges = assign_charges(atoms)
        hydrophobic = identify_hydrophobic_residues(chain_residues)
        
        if incremental_save:
            with open(f"{output_prefix}_charges.pkl", 'wb') as f:
                pickle.dump(charges, f)
            
            with open(f"{output_prefix}_hydrophobic.pkl", 'wb') as f:
                pickle.dump(hydrophobic, f)
            
            # Create index file
            index = {
                "atoms": f"{output_prefix}_atoms.pkl",
                "backbone": f"{output_prefix}_backbone.pkl",
                "secondary_structure": f"{output_prefix}_ss.pkl",
                "edge_pairs": f"{output_prefix}_edge_pairs.pkl",
                "bond_lengths": f"{output_prefix}_bond_lengths.pkl",
                "angles": f"{output_prefix}_angles.pkl",
                "torsions": f"{output_prefix}_torsions.pkl",
                "charges": f"{output_prefix}_charges.pkl",
                "hydrophobic": f"{output_prefix}_hydrophobic.pkl"
            }
            
            with open(f"{output_prefix}_index.json", 'w') as f:
                json.dump(index, f)
            
            return None
        
        # If not incremental save, return all data
        return {
            'backbone_atoms': backbone,
            'secondary_structure': chain_ss,
            'edge_indices': edge_index_pairs,
            'bond_lengths': bond_lengths,
            'angles': angles,
            'torsions': torsions,
            'charges': charges,
            'hydrophobic_residues': hydrophobic
        }
            
    except Exception as e:
        print(f"  Error processing chain {chain_id}: {str(e)}")
        traceback.print_exc()
        return {} if not incremental_save else None

def load_chain_data(chain_output_prefix):
    """
    Load chain data from incremental files.
    
    Parameters:
    -----------
    chain_output_prefix : str
        Prefix used for chain files
        
    Returns:
    --------
    dict
        Combined chain data
    """
    chain_data = {}
    
    try:
        # Try to load index file
        with open(f"{chain_output_prefix}_index.json", 'r') as f:
            index = json.load(f)
        
        # Load each file
        if 'backbone' in index:
            with open(index['backbone'], 'rb') as f:
                chain_data['backbone_atoms'] = pickle.load(f)
        
        if 'secondary_structure' in index:
            with open(index['secondary_structure'], 'rb') as f:
                chain_data['secondary_structure'] = pickle.load(f)
        
        if 'edge_pairs' in index:
            with open(index['edge_pairs'], 'rb') as f:
                chain_data['edge_indices'] = pickle.load(f)
        
        if 'bond_lengths' in index:
            with open(index['bond_lengths'], 'rb') as f:
                chain_data['bond_lengths'] = pickle.load(f)
        
        if 'angles' in index:
            with open(index['angles'], 'rb') as f:
                chain_data['angles'] = pickle.load(f)
        
        if 'torsions' in index:
            with open(index['torsions'], 'rb') as f:
                chain_data['torsions'] = pickle.load(f)
        
        if 'charges' in index:
            with open(index['charges'], 'rb') as f:
                chain_data['charges'] = pickle.load(f)
        
        if 'hydrophobic' in index:
            with open(index['hydrophobic'], 'rb') as f:
                chain_data['hydrophobic_residues'] = pickle.load(f)
        
        return chain_data
        
    except Exception as e:
        print(f"  Error loading chain data: {str(e)}")
        return {}

# Helper functions optimized for memory efficiency

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

def extract_backbone_atoms(chain_atoms):
    """Extract backbone atoms (N, CA, C, O) from a chain."""
    backbone = {}
    
    for atom in chain_atoms:
        chain_id, res_id, res_name, atom_name, element, position = atom
        
        if atom_name in ['N', 'CA', 'C', 'O']:
            if res_id not in backbone:
                backbone[res_id] = {}
            backbone[res_id][atom_name] = position
    
    return backbone

def create_edge_index_memory_efficient(chain_atoms, positions, mode='pairs', distance_cutoff=5.0):
    """
    Create edge index with memory efficiency in mind.
    Uses distance-based cutoff to reduce the number of edges.
    
    Parameters:
    -----------
    chain_atoms : list
        List of atom data
    positions : ndarray
        Array of atom positions
    mode : str
        'pairs', 'triplets', or 'quadruplets'
    distance_cutoff : float
        Maximum distance between atoms to create an edge
        
    Returns:
    --------
    np.ndarray or dict
        Edge index array or dictionary
    """
    n_atoms = len(chain_atoms)
    
    if mode == 'pairs':
        edges = []
        
        # Process in batches
        batch_size = 1000
        
        for i in range(0, n_atoms, batch_size):
            i_end = min(i + batch_size, n_atoms)
            
            for i_batch in range(i, i_end):
                res_i = chain_atoms[i_batch][1]  # Residue ID
                pos_i = positions[i_batch]
                
                for j in range(i_batch+1, n_atoms):
                    res_j = chain_atoms[j][1]  # Residue ID
                    pos_j = positions[j]
                    
                    # Check residue-based connectivity (same or adjacent)
                    residue_connected = (res_i == res_j or abs(res_i - res_j) == 1)
                    
                    # Check distance
                    if residue_connected:
                        # For residue-connected atoms, always add edge
                        edges.append([i_batch, j])
                        edges.append([j, i_batch])  # For undirected graph
                    else:
                        # For distant residues, only add if atoms are close
                        distance = np.linalg.norm(pos_j - pos_i)
                        if distance <= distance_cutoff:
                            edges.append([i_batch, j])
                            edges.append([j, i_batch])
            
            # Force garbage collection after each batch
            gc.collect()
        
        return np.array(edges).T if edges else np.zeros((2, 0), dtype=int)
    
    else:
        raise ValueError(f"Mode {mode} not implemented in memory-efficient version")

def calculate_bond_lengths_efficient(chain_atoms, edge_index):
    """Calculate bond lengths with memory efficiency."""
    bond_lengths = {}
    
    if isinstance(edge_index, np.ndarray) and edge_index.size > 0:
        # Process in batches
        batch_size = 5000
        total_edges = edge_index.shape[1]
        
        for start_idx in range(0, total_edges, batch_size):
            end_idx = min(start_idx + batch_size, total_edges)
            
            for idx in range(start_idx, end_idx):
                i, j = edge_index[0, idx], edge_index[1, idx]
                
                # Skip if already calculated (undirected graph)
                if (j, i) in bond_lengths:
                    bond_lengths[(i, j)] = bond_lengths[(j, i)]
                    continue
                
                # Calculate bond length
                pos_i = np.array(chain_atoms[i][5])
                pos_j = np.array(chain_atoms[j][5])
                bond_length = np.linalg.norm(pos_i - pos_j)
                
                bond_lengths[(i, j)] = bond_length
            
            # Force garbage collection after each batch
            gc.collect()
    
    return bond_lengths


def calculate_angles_memory_efficient(chain_atoms, positions, distance_cutoff=5.0, backbone_only=False):
    """
    Calculate angles with memory efficiency.
    
    This version doesn't create the full edge index first,
    but directly identifies valid triplets during calculation.
    
    Parameters:
        chain_atoms: List of atom information (chain_id, res_id, res_name, atom_name, element, position)
        positions: Numpy array of atom coordinates
        distance_cutoff: Maximum distance between atoms to be considered connected
        backbone_only: If True, only consider backbone atoms (N, CA, C, O)
    """
    angles = {}
    n_atoms = len(chain_atoms)
    
    if backbone_only:
        # Approach 1: Focus on backbone angles only
        # Group atoms by residue
        residue_groups = {}
        
        for i, atom in enumerate(chain_atoms):
            chain_id, res_id, res_name, atom_name, element, position = atom
            
            if atom_name in ['N', 'CA', 'C', 'O']:
                if res_id not in residue_groups:
                    residue_groups[res_id] = {}
                residue_groups[res_id][atom_name] = i
        
        # Sort residues by ID
        sorted_res_ids = sorted(residue_groups.keys())
        
        # Calculate angles within each residue
        triplet_count = 0
        
        # Intra-residue angles
        for res_id in sorted_res_ids:
            atoms = residue_groups[res_id]
            
            # Define typical backbone angles to calculate
            angle_patterns = [
                ('N', 'CA', 'C'),    # N-CA-C
                ('CA', 'C', 'O'),    # CA-C-O
                ('O', 'C', 'N'),     # O-C-N (if next residue exists)
                ('C', 'N', 'CA')     # C-N-CA (if previous residue exists)
            ]
            
            for i_name, j_name, k_name in angle_patterns:
                if i_name in atoms and j_name in atoms and k_name in atoms:
                    i = atoms[i_name]
                    j = atoms[j_name]
                    k = atoms[k_name]
                    
                    # Calculate angle
                    pos_i = positions[i]
                    pos_j = positions[j]
                    pos_k = positions[k]
                    
                    vec_ij = pos_i - pos_j
                    vec_jk = pos_k - pos_j
                    
                    # Compute angle
                    norm_ij = np.linalg.norm(vec_ij)
                    norm_jk = np.linalg.norm(vec_jk)
                    
                    if norm_ij < 1e-6 or norm_jk < 1e-6:
                        continue  # Skip if vectors too small
                    
                    cos_angle = np.dot(vec_ij, vec_jk) / (norm_ij * norm_jk)
                    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
                    
                    angles[(i, j, k)] = np.degrees(angle)
                    triplet_count += 1
        
        # Inter-residue angles
        for i in range(len(sorted_res_ids) - 1):
            curr_res = sorted_res_ids[i]
            next_res = sorted_res_ids[i + 1]
            
            # Skip if residues aren't sequential
            if next_res - curr_res != 1:
                continue
                
            curr_atoms = residue_groups[curr_res]
            next_atoms = residue_groups[next_res]
            
            # C-N-CA angle (connecting two residues)
            if 'C' in curr_atoms and 'N' in next_atoms and 'CA' in next_atoms:
                i = curr_atoms['C']
                j = next_atoms['N']
                k = next_atoms['CA']
                
                # Calculate angle
                pos_i = positions[i]
                pos_j = positions[j]
                pos_k = positions[k]
                
                vec_ij = pos_i - pos_j
                vec_jk = pos_k - pos_j
                
                norm_ij = np.linalg.norm(vec_ij)
                norm_jk = np.linalg.norm(vec_jk)
                
                if norm_ij >= 1e-6 and norm_jk >= 1e-6:
                    cos_angle = np.dot(vec_ij, vec_jk) / (norm_ij * norm_jk)
                    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
                    
                    angles[(i, j, k)] = np.degrees(angle)
                    triplet_count += 1
                    
            # CA-C-N angle (connecting two residues)
            if 'CA' in curr_atoms and 'C' in curr_atoms and 'N' in next_atoms:
                i = curr_atoms['CA']
                j = curr_atoms['C']
                k = next_atoms['N']
                
                # Calculate angle
                pos_i = positions[i]
                pos_j = positions[j]
                pos_k = positions[k]
                
                vec_ij = pos_i - pos_j
                vec_jk = pos_k - pos_j
                
                norm_ij = np.linalg.norm(vec_ij)
                norm_jk = np.linalg.norm(vec_jk)
                
                if norm_ij >= 1e-6 and norm_jk >= 1e-6:
                    cos_angle = np.dot(vec_ij, vec_jk) / (norm_ij * norm_jk)
                    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
                    
                    angles[(i, j, k)] = np.degrees(angle)
                    triplet_count += 1
        
        print(f"    Calculated {triplet_count} backbone angles")
        return angles
    
    else:
        # Approach 2: Original function for all angles
        # Find connected atom pairs (i-j and j-k must be connected)
        connected_pairs = set()
        
        # First identify connected pairs (atoms within same residue or sequential residues)
        for i in range(n_atoms):
            res_i = chain_atoms[i][1]
            pos_i = positions[i]
            
            for j in range(i+1, n_atoms):
                res_j = chain_atoms[j][1]
                pos_j = positions[j]
                
                # Check if connected by residue
                residue_connected = (res_i == res_j or abs(res_i - res_j) == 1)
                
                # Check distance
                distance = np.linalg.norm(pos_j - pos_i)
                
                if residue_connected or distance <= distance_cutoff:
                    connected_pairs.add((i, j))
                    connected_pairs.add((j, i))  # Undirected
        
        # Process in batches
        batch_size = 1000
        triplet_count = 0
        
        for j in range(n_atoms):
            # Find all potential triplets i-j-k where j is the central atom
            i_candidates = []
            k_candidates = []
            
            for i in range(n_atoms):
                if i != j and (i, j) in connected_pairs:
                    i_candidates.append(i)
            
            for k in range(n_atoms):
                if k != j and (j, k) in connected_pairs:
                    k_candidates.append(k)
            
            # Process candidates in batches
            for i_start in range(0, len(i_candidates), batch_size):
                i_end = min(i_start + batch_size, len(i_candidates))
                
                for k_start in range(0, len(k_candidates), batch_size):
                    k_end = min(k_start + batch_size, len(k_candidates))
                    
                    for i_idx in range(i_start, i_end):
                        i = i_candidates[i_idx]
                        
                        for k_idx in range(k_start, k_end):
                            k = k_candidates[k_idx]
                            
                            # Ensure i ≠ k
                            if i == k:
                                continue
                            
                            # Calculate angle
                            pos_i = positions[i]
                            pos_j = positions[j]
                            pos_k = positions[k]
                            
                            vec_ij = pos_i - pos_j
                            vec_jk = pos_k - pos_j
                            
                            # Compute angle
                            norm_ij = np.linalg.norm(vec_ij)
                            norm_jk = np.linalg.norm(vec_jk)
                            
                            if norm_ij < 1e-6 or norm_jk < 1e-6:
                                continue  # Skip if vectors too small
                            
                            cos_angle = np.dot(vec_ij, vec_jk) / (norm_ij * norm_jk)
                            angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
                            
                            angles[(i, j, k)] = np.degrees(angle)
                            triplet_count += 1
                    
                    # Force garbage collection after each batch
                    gc.collect()
        
        print(f"    Calculated {triplet_count} angles")
        return angles


def calculate_torsions_memory_efficient(chain_atoms, positions, distance_cutoff=2.0, backbone_only=True):
    """
    Calculate torsions with memory efficiency.
    
    Parameters:
    -----------
    chain_atoms : list
        List of atom data
    positions : ndarray
        Array of atom positions
    distance_cutoff : float
        Maximum distance to consider for connections
    backbone_only : bool
        If True, only calculate torsions for backbone atoms
        
    Returns:
    --------
    dict
        Dictionary of torsion angles
    """
    torsions = {}
    n_atoms = len(chain_atoms)
    
    if backbone_only:
        # Approach 1: Focus on backbone torsions only
        # Group atoms by residue
        residue_groups = {}
        
        for i, atom in enumerate(chain_atoms):
            chain_id, res_id, res_name, atom_name, element, position = atom
            
            if atom_name in ['N', 'CA', 'C', 'O']:
                if res_id not in residue_groups:
                    residue_groups[res_id] = {}
                residue_groups[res_id][atom_name] = i
        
        # Sort residues by ID
        sorted_res_ids = sorted(residue_groups.keys())
        
        # Calculate phi, psi, omega angles between consecutive residues
        for i in range(1, len(sorted_res_ids)):
            prev_res = sorted_res_ids[i-1]
            curr_res = sorted_res_ids[i]
            
            # Skip if residues aren't sequential
            if curr_res - prev_res != 1:
                continue
                
            prev_atoms = residue_groups.get(prev_res, {})
            curr_atoms = residue_groups.get(curr_res, {})
            
            # Calculate phi (C-N-CA-C)
            if ('C' in prev_atoms and 'N' in curr_atoms and 
                'CA' in curr_atoms and 'C' in curr_atoms):
                i = prev_atoms['C']
                j = curr_atoms['N']
                k = curr_atoms['CA']
                l = curr_atoms['C']
                
                torsion = calculate_single_torsion(
                    positions[i], positions[j], positions[k], positions[l])
                if torsion is not None:
                    torsions[(i, j, k, l)] = torsion
            
            # Calculate psi if next residue exists
            if i < len(sorted_res_ids) - 1:
                next_res = sorted_res_ids[i+1]
                
                # Skip if residues aren't sequential
                if next_res - curr_res != 1:
                    continue
                    
                next_atoms = residue_groups.get(next_res, {})
                
                if ('N' in curr_atoms and 'CA' in curr_atoms and 
                    'C' in curr_atoms and 'N' in next_atoms):
                    i = curr_atoms['N']
                    j = curr_atoms['CA']
                    k = curr_atoms['C']
                    l = next_atoms['N']
                    
                    torsion = calculate_single_torsion(
                        positions[i], positions[j], positions[k], positions[l])
                    if torsion is not None:
                        torsions[(i, j, k, l)] = torsion
            
            # Calculate omega
            if ('CA' in prev_atoms and 'C' in prev_atoms and 
                'N' in curr_atoms and 'CA' in curr_atoms):
                i = prev_atoms['CA']
                j = prev_atoms['C']
                k = curr_atoms['N']
                l = curr_atoms['CA']
                
                torsion = calculate_single_torsion(
                    positions[i], positions[j], positions[k], positions[l])
                if torsion is not None:
                    torsions[(i, j, k, l)] = torsion
        
        print(f"    Calculated {len(torsions)} backbone torsions")
        return torsions
    
    else:
        # Approach 2: Original function for all torsions
        # First identify connected pairs (similar to angle calculation)
        connected_pairs = set()
        
        for i in range(n_atoms):
            res_i = chain_atoms[i][1]
            pos_i = positions[i]
            
            for j in range(i+1, n_atoms):
                res_j = chain_atoms[j][1]
                pos_j = positions[j]
                
                # Check if connected by residue
                residue_connected = (res_i == res_j or abs(res_i - res_j) == 1)
                
                # Check distance
                distance = np.linalg.norm(pos_j - pos_i)
                
                if residue_connected or distance <= distance_cutoff:
                    connected_pairs.add((i, j))
                    connected_pairs.add((j, i))  # Undirected
        
        # Find valid quadruplets (i-j-k-l where each adjacent pair is connected)
        batch_size = 500
        quadruplet_count = 0
        
        # First find all potential j-k pairs
        jk_pairs = []
        for j in range(n_atoms):
            for k in range(j+1, n_atoms):
                if (j, k) in connected_pairs:
                    jk_pairs.append((j, k))
        
        # Process j-k pairs in batches
        for jk_start in range(0, len(jk_pairs), batch_size):
            jk_end = min(jk_start + batch_size, len(jk_pairs))
            
            for jk_idx in range(jk_start, jk_end):
                j, k = jk_pairs[jk_idx]
                
                # Find i's connected to j
                i_candidates = []
                for i in range(n_atoms):
                    if i != j and i != k and (i, j) in connected_pairs:
                        i_candidates.append(i)
                
                # Find l's connected to k
                l_candidates = []
                for l in range(n_atoms):
                    if l != i and l != j and l != k and (k, l) in connected_pairs:
                        l_candidates.append(l)
                
                # Process i-l candidates in smaller batches
                for i_start in range(0, len(i_candidates), batch_size):
                    i_end = min(i_start + batch_size, len(i_candidates))
                    
                    for l_start in range(0, len(l_candidates), batch_size):
                        l_end = min(l_start + batch_size, len(l_candidates))
                        
                        for i_idx in range(i_start, i_end):
                            i = i_candidates[i_idx]
                            
                            for l_idx in range(l_start, l_end):
                                l = l_candidates[l_idx]
                                
                                # Calculate torsion
                                torsion = calculate_single_torsion(
                                    positions[i], positions[j], positions[k], positions[l])
                                if torsion is not None:
                                    torsions[(i, j, k, l)] = torsion
                                    quadruplet_count += 1
                        
                        # Force garbage collection
                        gc.collect()
        
        print(f"    Calculated {quadruplet_count} torsions")
        return torsions

def calculate_single_torsion(pos_i, pos_j, pos_k, pos_l):
    """Calculate a single torsion angle from 4 positions."""
    # Define vectors
    v1 = pos_j - pos_i
    v2 = pos_k - pos_j
    v3 = pos_l - pos_k
    
    # Skip if any vector is too short
    if (np.linalg.norm(v1) < 1e-6 or 
        np.linalg.norm(v2) < 1e-6 or 
        np.linalg.norm(v3) < 1e-6):
        return None
    
    # Cross products
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    
    # Skip if either normal is too short
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    
    if n1_norm < 1e-6 or n2_norm < 1e-6:
        return None
    
    # Normalize
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    
    # Compute torsion angle
    v2_norm = np.linalg.norm(v2)
    torsion = np.arctan2(
        np.dot(np.cross(n1, n2), v2 / v2_norm), 
        np.dot(n1, n2)
    )
    
    return np.degrees(torsion)

def calculate_torsions_all_atoms(chain_atoms, positions, distance_cutoff=2.0):
    """
    Calculate torsions with memory efficiency.
    
    This version doesn't create the full edge index first,
    but directly identifies valid quadruplets during calculation.
    """
    torsions = {}
    n_atoms = len(chain_atoms)
    
    # First identify connected pairs (similar to angle calculation)
    connected_pairs = set()
    
    for i in range(n_atoms):
        res_i = chain_atoms[i][1]
        pos_i = positions[i]
        
        for j in range(i+1, n_atoms):
            res_j = chain_atoms[j][1]
            pos_j = positions[j]
            
            # Check if connected by residue
            residue_connected = (res_i == res_j or abs(res_i - res_j) == 1)
            
            # Check distance
            distance = np.linalg.norm(pos_j - pos_i)
            
            if residue_connected or distance <= distance_cutoff:
                connected_pairs.add((i, j))
                connected_pairs.add((j, i))  # Undirected
    
    # Find valid quadruplets (i-j-k-l where each adjacent pair is connected)
    batch_size = 500
    quadruplet_count = 0
    
    # First find all potential j-k pairs
    jk_pairs = []
    for j in range(n_atoms):
        for k in range(j+1, n_atoms):
            if (j, k) in connected_pairs:
                jk_pairs.append((j, k))
    
    # Process j-k pairs in batches
    for jk_start in range(0, len(jk_pairs), batch_size):
        jk_end = min(jk_start + batch_size, len(jk_pairs))
        
        for jk_idx in range(jk_start, jk_end):
            j, k = jk_pairs[jk_idx]
            
            # Find i's connected to j
            i_candidates = []
            for i in range(n_atoms):
                if i != j and i != k and (i, j) in connected_pairs:
                    i_candidates.append(i)
            
            # Find l's connected to k
            l_candidates = []
            for l in range(n_atoms):
                if l != i and l != j and l != k and (k, l) in connected_pairs:
                    l_candidates.append(l)
            
            # Process i-l candidates in smaller batches
            for i_start in range(0, len(i_candidates), batch_size):
                i_end = min(i_start + batch_size, len(i_candidates))
                
                for l_start in range(0, len(l_candidates), batch_size):
                    l_end = min(l_start + batch_size, len(l_candidates))
                    
                    for i_idx in range(i_start, i_end):
                        i = i_candidates[i_idx]
                        
                        for l_idx in range(l_start, l_end):
                            l = l_candidates[l_idx]
                            
                            # Calculate torsion
                            pos_i = positions[i]
                            pos_j = positions[j]
                            pos_k = positions[k]
                            pos_l = positions[l]
                            
                            # Define vectors
                            v1 = pos_j - pos_i
                            v2 = pos_k - pos_j
                            v3 = pos_l - pos_k
                            
                            # Skip if any vector is too short
                            if (np.linalg.norm(v1) < 1e-6 or 
                                np.linalg.norm(v2) < 1e-6 or 
                                np.linalg.norm(v3) < 1e-6):
                                continue
                            
                            # Cross products
                            n1 = np.cross(v1, v2)
                            n2 = np.cross(v2, v3)
                            
                            # Skip if either normal is too short
                            n1_norm = np.linalg.norm(n1)
                            n2_norm = np.linalg.norm(n2)
                            
                            if n1_norm < 1e-6 or n2_norm < 1e-6:
                                continue
                            
                            # Normalize
                            n1 = n1 / n1_norm
                            n2 = n2 / n2_norm
                            
                            # Compute torsion angle
                            v2_norm = np.linalg.norm(v2)
                            torsion = np.arctan2(
                                np.dot(np.cross(n1, n2), v2 / v2_norm), 
                                np.dot(n1, n2)
                            )
                            
                            torsions[(i, j, k, l)] = np.degrees(torsion)
                            quadruplet_count += 1
                    
                    # Force garbage collection
                    gc.collect()
    
    print(f"    Calculated {quadruplet_count} torsions")
    return torsions

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

def create_protein_graph(pdb_path, ss_data=None):
    """
    Create a protein graph using Graphein or fallback to custom implementation.
    
    Parameters:
    pdb_path -- Path to the PDB file
    ss_data -- Optional secondary structure data dictionary
    
    Returns:
    NetworkX graph
    """
    try:
        # Try the current Graphein API
        import graphein.protein as gp
        
        # Check the actual function signature
        import inspect
        sig = inspect.signature(gp.construct_graph)
        
        # Check if 'pdb_path' is a parameter
        if 'pdb_path' in sig.parameters:
            graph = gp.construct_graph(pdb_path=pdb_path)
        elif 'path' in sig.parameters:
            # Try with 'path' instead
            graph = gp.construct_graph(path=pdb_path)
        else:
            # Check if the function expects the first argument to be the path
            params = list(sig.parameters.keys())
            if len(params) > 0:
                # Try with positional argument
                graph = gp.construct_graph(pdb_path)
            else:
                raise ValueError("Cannot determine correct parameters for construct_graph")
                
        return graph
        
    except Exception as e:
        print(f"Graphein error: {str(e)}")
        print("Falling back to custom graph implementation...")
        
        # Create a simple graph using NetworkX and BioPython
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        
        # Create basic graph
        graph = nx.Graph()
        
        # Add nodes for each residue
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    # Skip non-standard residues
                    if residue.get_id()[0] != " ":
                        continue
                        
                    res_id = residue.get_id()[1]
                    res_name = residue.get_resname()
                    
                    # Add node
                    node_id = f"{chain_id}:{res_name}:{res_id}"
                    graph.add_node(
                        node_id,
                        chain_id=chain_id,
                        residue_number=res_id,
                        residue_name=res_name
                    )
                    
                    # Add atom data and track if backbone atoms are present
                    has_backbone = False
                    for atom in residue:
                        atom_name = atom.get_name()
                        coord = atom.get_coord()
                        
                        # Store coordinates of backbone atoms
                        if atom_name in ['CA', 'N', 'C', 'O']:
                            graph.nodes[node_id][atom_name] = coord
                            has_backbone = True
                    
                    # Set backbone flag
                    graph.nodes[node_id]['has_backbone'] = has_backbone
                    
                    # Store CA coordinates for distance calculations
                    if 'CA' in graph.nodes[node_id]:
                        graph.nodes[node_id]['coords'] = graph.nodes[node_id]['CA']
                    
                    # Add secondary structure if available
                    if ss_data and (chain_id, res_id) in ss_data:
                        graph.nodes[node_id]['ss'] = ss_data[(chain_id, res_id)]
                    else:
                        graph.nodes[node_id]['ss'] = '?'
        
        # Add edges between consecutive residues and close residues
        # First, sort nodes by chain and residue number
        sorted_nodes = sorted(graph.nodes(), key=lambda x: (x.split(':')[0], int(x.split(':')[2]) if x.split(':')[2].isdigit() else 0))
        
        # Add backbone connections
        for i in range(len(sorted_nodes) - 1):
            node1 = sorted_nodes[i]
            chain1, _, res1 = node1.split(':')
            
            for j in range(i + 1, min(i + 20, len(sorted_nodes))):  # Check only nearby residues
                node2 = sorted_nodes[j]
                chain2, _, res2 = node2.split(':')
                
                # Add edge if same chain and sequential
                if chain1 == chain2 and abs(int(res1) - int(res2)) == 1:
                    graph.add_edge(node1, node2, kind={'peptide_bond'})
                
                # Add edge if CA atoms are within 7 Angstroms
                if 'coords' in graph.nodes[node1] and 'coords' in graph.nodes[node2]:
                    ca1 = graph.nodes[node1]['coords']
                    ca2 = graph.nodes[node2]['coords']
                    distance = np.linalg.norm(ca1 - ca2)
                    
                    if distance < 7.0:  # 7 Angstrom cutoff
                        graph.add_edge(node1, node2, kind={'contact'}, distance=distance)
        
        return graph


def add_structure_to_graph(graph, ss_data):
    """
    Add secondary structure and backbone information to a protein graph.
    
    Parameters:
    graph -- NetworkX graph to update
    ss_data -- Dictionary of secondary structure data
    
    Returns:
    Updated NetworkX graph
    """
    for node, data in graph.nodes(data=True):
        # Extract chain and residue information
        if ':' in node:  # Format like 'H:VAL:1'
            parts = node.split(':')
            chain_id = parts[0]
            if len(parts) >= 3:
                try:
                    residue_id = int(parts[2])
                except ValueError:
                    residue_id = parts[2]
        else:
            # Try to extract from attributes
            chain_id = data.get('chain_id', '')
            residue_id = data.get('residue_number', None)
        
        # Add secondary structure if available
        if (chain_id, residue_id) in ss_data:
            graph.nodes[node]['ss'] = ss_data[(chain_id, residue_id)]
        else:
            graph.nodes[node]['ss'] = '?'
        
        # Check for backbone atoms using exact attribute names from Graphein
        backbone_keys = ['CA', 'C', 'N', 'O', 'atom_type', 'element_symbol', 'coords']
        backbone_atoms_present = False
        
        for key in data.keys():
            # Check for backbone atom coordinates or any CA atom indicator
            if key in ['CA', 'C', 'N', 'O']:
                backbone_atoms_present = True
                break
            elif key == 'atom_type' and data[key] in ['CA', 'C', 'N', 'O']:
                backbone_atoms_present = True
                break
        
        # For debugging, print the first node's data
        if node == list(graph.nodes())[0]:
            print(f"First node data keys: {list(data.keys())}")
            print(f"Checking for backbone atoms: {backbone_atoms_present}")
            
        # Add backbone information
        graph.nodes[node]['has_backbone'] = backbone_atoms_present
    
    return graph

def add_backbone_from_pdb(graph, pdb_path):
    """
    Add backbone information to the graph by directly extracting it from the PDB file.
    
    Parameters:
    graph -- NetworkX graph to update
    pdb_path -- Path to the PDB file
    
    Returns:
    Updated NetworkX graph
    """
    # Parse PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    # Extract backbone atoms
    backbone_info = {}
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                if residue.get_id()[0] != " ":  # Skip non-standard residues
                    continue
                    
                res_id = residue.get_id()[1]
                
                # Check for backbone atoms
                backbone_atoms = {}
                for atom_name in ['CA', 'N', 'C', 'O']:
                    if atom_name in residue:
                        backbone_atoms[atom_name] = residue[atom_name].get_coord()
                
                if backbone_atoms:
                    backbone_info[(chain_id, res_id)] = backbone_atoms
    
    # Add backbone info to graph
    for node, data in graph.nodes(data=True):
        # Extract chain and residue ID
        if ':' in node:
            parts = node.split(':')
            chain_id = parts[0]
            if len(parts) >= 3:
                try:
                    res_id = int(parts[2])
                except ValueError:
                    res_id = parts[2]
                    
                # Check if we have backbone info for this residue
                if (chain_id, res_id) in backbone_info:
                    # Add backbone atom coordinates
                    for atom_name, coord in backbone_info[(chain_id, res_id)].items():
                        graph.nodes[node][atom_name] = coord
                    
                    # Set backbone flag
                    graph.nodes[node]['has_backbone'] = True
                else:
                    graph.nodes[node]['has_backbone'] = False
        else:
            # Extract from attributes if possible
            chain_id = data.get('chain_id', '')
            res_id = data.get('residue_number', None)
            
            if chain_id and res_id and (chain_id, res_id) in backbone_info:
                # Add backbone atom coordinates
                for atom_name, coord in backbone_info[(chain_id, res_id)].items():
                    graph.nodes[node][atom_name] = coord
                
                # Set backbone flag
                graph.nodes[node]['has_backbone'] = True
    
    return graph

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












def visualize_protein_graph(graph_path, output_path=None, show=True, color_by='ss', node_size=100, figsize=(12, 10)):
    """
    Visualize a protein graph from a saved NetworkX pickle file.
    
    Parameters:
    -----------
    graph_path : str
        Path to the saved NetworkX graph pickle file
    output_path : str, optional
        Path to save the visualization image
    show : bool
        Whether to display the graph
    color_by : str
        Property to color nodes by ('ss' for secondary structure,
        'chain', 'residue_name', or any other node attribute)
    node_size : int
        Size of nodes in the visualization
    figsize : tuple
        Figure size (width, height) in inches
    
    Returns:
    --------
    matplotlib figure
    """
    import pickle
    import networkx as nx
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Load the graph
    with open(graph_path, 'rb') as f:
        G = pickle.load(f)
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Create position layout
    # Use 3D positions if available, otherwise spring layout
    pos = {}
    for node, data in G.nodes(data=True):
        if 'coords' in data:
            # Use first two coordinates for 2D visualization
            pos[node] = data['coords'][:2]
        elif 'CA' in data:
            # If explicit CA coordinates exist
            pos[node] = data['CA'][:2]
    
    # If we couldn't get positions from the graph, use spring layout
    if not pos:
        pos = nx.spring_layout(G, seed=42)
    
    # Color mapping based on the selected property
    if color_by == 'ss':
        # Secondary structure coloring
        ss_types = {'H': 'red', 'E': 'yellow', 'B': 'green', 
                    'G': 'orange', 'I': 'purple', 'T': 'cyan', 
                    'S': 'blue', '?': 'gray', ' ': 'gray'}
        
        node_colors = []
        for node, data in G.nodes(data=True):
            ss = data.get('ss', '?')
            node_colors.append(ss_types.get(ss, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for secondary structure
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=ss) 
                          for ss, color in ss_types.items() if any(data.get('ss', '?') == ss for _, data in G.nodes(data=True))]
        plt.legend(handles=legend_elements, title="Secondary Structure")
        
    elif color_by == 'chain':
        # Color by chain ID
        chains = set()
        for node in G.nodes():
            if ':' in node:
                chains.add(node.split(':')[0])
            else:
                chains.add(G.nodes[node].get('chain_id', 'Unknown'))
        
        # Create color map
        cmap = plt.cm.get_cmap('tab10', len(chains))
        chain_colors = {chain: cmap(i) for i, chain in enumerate(sorted(chains))}
        
        node_colors = []
        for node in G.nodes():
            if ':' in node:
                chain = node.split(':')[0]
            else:
                chain = G.nodes[node].get('chain_id', 'Unknown')
            node_colors.append(chain_colors.get(chain, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for chains
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=chain) 
                          for chain, color in chain_colors.items()]
        plt.legend(handles=legend_elements, title="Chains")
        
    elif color_by == 'residue_name':
        # Color by amino acid type
        aa_types = set()
        for node, data in G.nodes(data=True):
            if ':' in node:
                aa_types.add(node.split(':')[1])
            else:
                aa_types.add(data.get('residue_name', 'Unknown'))
        
        # Create color map
        cmap = plt.cm.get_cmap('tab20', len(aa_types))
        aa_colors = {aa: cmap(i) for i, aa in enumerate(sorted(aa_types))}
        
        node_colors = []
        for node in G.nodes():
            if ':' in node:
                aa = node.split(':')[1]
            else:
                aa = G.nodes[node].get('residue_name', 'Unknown')
            node_colors.append(aa_colors.get(aa, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for amino acids
        from matplotlib.patches import Patch
        # Only show most common AAs in legend to avoid overcrowding
        aa_counts = {}
        for node in G.nodes():
            if ':' in node:
                aa = node.split(':')[1]
            else:
                aa = G.nodes[node].get('residue_name', 'Unknown')
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        # Get top amino acids
        top_aas = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        legend_elements = [Patch(facecolor=aa_colors[aa], label=f"{aa} ({count})") 
                          for aa, count in top_aas]
        plt.legend(handles=legend_elements, title="Top Amino Acids")
    
    else:
        # Try to color by any other attribute
        try:
            values = [G.nodes[node].get(color_by, 0) for node in G.nodes()]
            
            if all(isinstance(v, (int, float)) for v in values if v is not None):
                # Numeric values - use colormap
                cmap = plt.cm.viridis
                vmin = min(v for v in values if v is not None)
                vmax = max(v for v in values if v is not None)
                
                # Replace None with vmin
                values = [vmin if v is None else v for v in values]
                
                nx.draw_networkx_nodes(G, pos, node_color=values, 
                                      cmap=cmap, vmin=vmin, vmax=vmax, node_size=node_size)
                plt.colorbar(label=color_by)
            else:
                # Categorical values
                unique_values = set(v for v in values if v is not None)
                cmap = plt.cm.get_cmap('tab20', len(unique_values))
                value_colors = {val: cmap(i) for i, val in enumerate(sorted(unique_values))}
                
                node_colors = [value_colors.get(G.nodes[node].get(color_by, None), 'gray') 
                              for node in G.nodes()]
                nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
                
                # Create legend
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor=color, label=str(val)) 
                                  for val, color in value_colors.items()]
                plt.legend(handles=legend_elements, title=color_by)
        except:
            # Fallback to default coloring
            nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=node_size)
    
    # Draw edges with alpha for clarity
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5)
    
    # Optional: Draw node labels for smaller graphs (uncomment if needed)
    # if G.number_of_nodes() < 100:  # Only show labels for smaller graphs
    #     nx.draw_networkx_labels(G, pos, font_size=8)
    
    # Add graph info
    plt.title(f"Protein Graph: {graph_path.split('/')[-1]}")
    plt.text(0.02, 0.02, f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}", 
             transform=plt.gca().transAxes, fontsize=10)
    
    plt.axis('off')
    
    # Save if requested
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
    
    # Show if requested
    if show:
        plt.show()
    else:
        plt.close()
    
    return plt.gcf()

def visualize_all_graphs(output_dir, visualization_dir=None):
    """
    Visualize all graphs in the output directory.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing processed PDB data
    visualization_dir : str, optional
        Directory to save visualizations
    """
    import os
    import glob
    
    if visualization_dir is None:
        visualization_dir = os.path.join(output_dir, "visualizations")
    
    os.makedirs(visualization_dir, exist_ok=True)
    
    # Find all graph files
    graph_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith("_graph.pkl"):
                graph_files.append(os.path.join(root, file))
    
    print(f"Found {len(graph_files)} graph files to visualize")
    
    for graph_file in graph_files:
        pdb_id = os.path.basename(graph_file).replace("_graph.pkl", "")
        print(f"Visualizing {pdb_id}...")
        
        # Create visualizations with different colorings
        for color_scheme in ['ss', 'chain', 'residue_name']:
            output_path = os.path.join(visualization_dir, f"{pdb_id}_{color_scheme}.png")
            
            visualize_protein_graph(
                graph_file,
                output_path=output_path,
                color_by=color_scheme,
                node_size=50,
                show=False  # Don't display, just save
            )
    
    print(f"All visualizations saved to {visualization_dir}")





if __name__ == "__main__":
    folder_path = "new_protein_data"  # Your PDB folder
    output_path = "new_protein_data_networkx_new"
    
    # Process with memory optimizations but full features
    summary = process_pdb_full_features_memory_efficient(
        folder_path,
        output_path=output_path,
        max_file_size_mb=25  # Skip files larger than this
    )


    # # After processing is complete
    # graph_path = "fabs_networkx/7JKB/7JKB_graph.pkl"
    # visualize_protein_graph(
    #     graph_path,
    #     output_path="7JKB_graph_visualization.png",
    #     color_by='ss',  # Color by secondary structure
    #     node_size=50    # Smaller nodes for better visibility
    # )

    # # Try other visualizations:
    # # Color by chain
    # visualize_protein_graph(
    #     graph_path,
    #     output_path="7JKB_graph_by_chain.png",
    #     color_by='chain',
    #     node_size=50
    # )

    # # Color by amino acid type
    # visualize_protein_graph(
    #     graph_path,
    #     output_path="7JKB_graph_by_residue.png",
    #     color_by='residue_name',
    #     node_size=50
    # )


    # visualize_all_graphs("fabs_networkx")