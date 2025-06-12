import numpy as np
import gc


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
