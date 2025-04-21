import os
import numpy as np
import networkx as nx
import pickle
import json
from tqdm import tqdm
import random
import time

def generate_synthetic_protein_graph(num_nodes=100, seed=None):
    """
    Generate a synthetic protein graph with specified number of nodes.

    Parameters:
    -----------
    num_nodes : int
        Number of nodes in the graph
    seed : int or None
        Random seed for reproducibility

    Returns:
    --------
    nx.Graph
        NetworkX graph representing a protein structure
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Create an empty graph
    G = nx.Graph()

    # Define chains (typically 2-3 chains in a protein)
    chains = ['A', 'B']

    # Assign residues to chains
    chain_sizes = []
    remaining_nodes = num_nodes

    for i in range(len(chains) - 1):
        # Random proportion of nodes to this chain
        size = int(remaining_nodes * random.uniform(0.3, 0.7))
        chain_sizes.append(size)
        remaining_nodes -= size

    chain_sizes.append(remaining_nodes)

    # 20 standard amino acids
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                   'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                   'THR', 'TRP', 'TYR', 'VAL']

    # Secondary structure types
    ss_types = ['H', 'E', 'B', 'G', 'T', 'S', '?']
    ss_weights = [0.4, 0.3, 0.05, 0.05, 0.1, 0.05, 0.05]  # More helices and sheets

    # Generate nodes for each chain
    current_node = 0
    nodes_by_chain = {}

    for i, chain_id in enumerate(chains):
        chain_size = chain_sizes[i]
        nodes_by_chain[chain_id] = []

        for j in range(chain_size):
            res_number = j + 1
            res_name = random.choice(amino_acids)

            # Create node ID
            node_id = f"{chain_id}:{res_name}:{res_number}"
            nodes_by_chain[chain_id].append(node_id)

            # Generate random 3D coordinates
            coords = np.random.normal(0, 10, 3)

            # Random secondary structure
            ss = random.choices(ss_types, weights=ss_weights)[0]

            # Add node with attributes
            G.add_node(
                node_id,
                chain_id=chain_id,
                residue_number=res_number,
                residue_name=res_name,
                ss=ss,
                has_backbone=True,
                coords=coords,
                # Backbone atom coordinates
                CA=coords,
                N=coords + np.random.normal(0, 1, 3),
                C=coords + np.random.normal(0, 1, 3),
                O=coords + np.random.normal(0, 1, 3)
            )

            current_node += 1

    # Add sequential edges (peptide bonds)
    for chain_id, node_list in nodes_by_chain.items():
        for i in range(len(node_list) - 1):
            G.add_edge(node_list[i], node_list[i+1], kind={'peptide_bond'}, distance=3.8)

    # Add spatial proximity edges
    for chain_id, node_list in nodes_by_chain.items():
        for i in range(len(node_list)):
            # Add some random long-range connections within the same chain
            for _ in range(random.randint(1, 3)):
                j = random.randint(0, len(node_list) - 1)
                if abs(i - j) > 3:  # Not too close in sequence
                    distance = np.random.uniform(4.0, 6.9)
                    G.add_edge(node_list[i], node_list[j], kind={'contact'}, distance=distance)

            # Add some inter-chain contacts
            for other_chain, other_nodes in nodes_by_chain.items():
                if other_chain != chain_id:
                    for _ in range(random.randint(0, 2)):
                        j = random.randint(0, len(other_nodes) - 1)
                        distance = np.random.uniform(4.0, 6.9)
                        G.add_edge(node_list[i], other_nodes[j], kind={'contact'}, distance=distance)

    return G

def generate_associated_data(graph, seed=None):
    """
    Generate associated data structures for a protein graph.

    Parameters:
    -----------
    graph : nx.Graph
        NetworkX graph representing a protein structure
    seed : int or None
        Random seed for reproducibility

    Returns:
    --------
    dict
        Dictionary containing associated data
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    # Extract all nodes, organize by chain
    residues_by_chain = {}
    for node, data in graph.nodes(data=True):
        chain_id = data['chain_id']
        res_id = data['residue_number']
        res_name = data['residue_name']

        if chain_id not in residues_by_chain:
            residues_by_chain[chain_id] = []

        residues_by_chain[chain_id].append((res_id, res_name))

    # Sort residues by ID within each chain
    for chain_id in residues_by_chain:
        residues_by_chain[chain_id].sort(key=lambda x: x[0])

    # Generate mock atoms list
    atoms = []
    atom_idx = 0
    atoms_by_chain = {}

    for chain_id, residues in residues_by_chain.items():
        atoms_by_chain[chain_id] = []

        for res_id, res_name in residues:
            # Basic backbone atoms
            for atom_name in ['N', 'CA', 'C', 'O']:
                element = atom_name[0]
                position = np.random.normal(0, 10, 3)
                atoms.append((chain_id, res_id, res_name, atom_name, element, position))
                atoms_by_chain[chain_id].append(atom_idx)
                atom_idx += 1

            # Some side chain atoms (simplified)
            num_sc_atoms = random.randint(1, 4)
            for i in range(num_sc_atoms):
                atom_name = f"CB{i+1}" if i > 0 else "CB"
                element = 'C'
                position = np.random.normal(0, 10, 3)
                atoms.append((chain_id, res_id, res_name, atom_name, element, position))
                atoms_by_chain[chain_id].append(atom_idx)
                atom_idx += 1

    # Generate other data structures
    data = {
        'pdb_id': 'synthetic',
        'atoms': atoms,
        'residues_by_chain': residues_by_chain,
        'relative_positions': {},
        'edge_indices': {},
        'bond_lengths': {},
        'angles': {},
        'torsions': {},
        'secondary_structure': {},
        'backbone_atoms': {},
        'charges': {},
        'hydrophobic_residues': {}
    }

    # Generate data for each chain
    for chain_id, chain_residues in residues_by_chain.items():
        # Extract chain atoms
        chain_atoms = [atom for atom in atoms if atom[0] == chain_id]
        chain_atom_indices = atoms_by_chain[chain_id]

        # Backbone atoms
        backbone = {}
        for res_id, _ in chain_residues:
            backbone[res_id] = {}
            for atom_name in ['N', 'CA', 'C', 'O']:
                backbone[res_id][atom_name] = np.random.normal(0, 10, 3)

        # Secondary structure (already in graph)
        ss = {}
        for node, node_data in graph.nodes(data=True):
            if ':' in node:
                parts = node.split(':')
                if parts[0] == chain_id:
                    try:
                        res_id = int(parts[2])
                        ss[res_id] = node_data.get('ss', '?')
                    except ValueError:
                        pass

        # Bond lengths (random but realistic)
        bond_lengths = {}
        num_bonds = len(chain_atom_indices) * 2  # Approx 2 bonds per atom
        for _ in range(num_bonds):
            i = random.choice(chain_atom_indices)
            j = random.choice(chain_atom_indices)
            if i != j and (j, i) not in bond_lengths:
                bond_lengths[(i, j)] = np.random.uniform(1.0, 2.0)  # Typical bond length range

        # Angles
        angles = {}
        num_angles = len(chain_atom_indices)  # Approx 1 angle per atom
        for _ in range(num_angles):
            i = random.choice(chain_atom_indices)
            j = random.choice(chain_atom_indices)
            k = random.choice(chain_atom_indices)
            if i != j and j != k and i != k and (k, j, i) not in angles:
                angles[(i, j, k)] = np.random.uniform(80, 150)  # Typical angle range

        # Torsions (focus on backbone)
        torsions = {}
        for res_idx in range(len(chain_residues) - 1):
            # Simplified: just create random torsions for consecutive residues
            # In real data, these would be phi/psi/omega angles
            torsions[(res_idx*4, res_idx*4+1, res_idx*4+2, res_idx*4+3)] = np.random.uniform(-180, 180)

        # Charges
        charges = {}
        for i in chain_atom_indices:
            charges[i] = np.random.uniform(-0.5, 0.5)

        # Hydrophobic residues
        hydrophobic_aas = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        hydrophobic = {res_id for res_id, res_name in chain_residues if res_name in hydrophobic_aas}

        # Positions and edge indices
        positions = np.array([atom[5] for atom in chain_atoms])
        # Create edge index pairs (simplified)
        edge_pairs = np.zeros((2, len(bond_lengths)), dtype=int)
        for idx, (i, j) in enumerate(bond_lengths.keys()):
            edge_pairs[0, idx] = i
            edge_pairs[1, idx] = j

        # Add all computed data to the main dictionary
        data['backbone_atoms'][chain_id] = backbone
        data['secondary_structure'][chain_id] = ss
        data['bond_lengths'][chain_id] = bond_lengths
        data['angles'][chain_id] = angles
        data['torsions'][chain_id] = torsions
        data['charges'][chain_id] = charges
        data['hydrophobic_residues'][chain_id] = hydrophobic
        data['edge_indices'][chain_id] = edge_pairs

    return data

def save_graph_and_data(graph, data, output_dir, pdb_id, save_individual_files=True):
    """
    Save graph and associated data to disk.

    Parameters:
    -----------
    graph : nx.Graph
        NetworkX graph
    data : dict
        Associated data dictionary
    output_dir : str
        Output directory
    pdb_id : str
        PDB identifier for filenaming
    save_individual_files : bool
        Whether to save individual chain files
    """
    # Create directories
    os.makedirs(output_dir, exist_ok=True)
    pdb_dir = os.path.join(output_dir, pdb_id)
    os.makedirs(pdb_dir, exist_ok=True)

    # Save graph
    graph_path = os.path.join(pdb_dir, f"{pdb_id}_graph.pkl")
    with open(graph_path, 'wb') as f:
        pickle.dump(graph, f)

    # Save main data dictionary
    data_path = os.path.join(pdb_dir, f"{pdb_id}_data.pkl")
    with open(data_path, 'wb') as f:
        pickle.dump(data, f)

    if save_individual_files:
        # Save individual chain files
        for chain_id in data['residues_by_chain'].keys():
            chain_dir = os.path.join(pdb_dir, f"{chain_id}")
            os.makedirs(chain_dir, exist_ok=True)

            # Chain atoms
            chain_atoms = [atom for atom in data['atoms'] if atom[0] == chain_id]
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_atoms.pkl"), 'wb') as f:
                pickle.dump(chain_atoms, f)

            # Backbone
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_backbone.pkl"), 'wb') as f:
                pickle.dump(data['backbone_atoms'].get(chain_id, {}), f)

            # Secondary structure
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_ss.pkl"), 'wb') as f:
                pickle.dump(data['secondary_structure'].get(chain_id, {}), f)

            # Bond lengths
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_bond_lengths.pkl"), 'wb') as f:
                pickle.dump(data['bond_lengths'].get(chain_id, {}), f)

            # Angles
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_angles.pkl"), 'wb') as f:
                pickle.dump(data['angles'].get(chain_id, {}), f)

            # Torsions
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_torsions.pkl"), 'wb') as f:
                pickle.dump(data['torsions'].get(chain_id, {}), f)

            # Charges
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_charges.pkl"), 'wb') as f:
                pickle.dump(data['charges'].get(chain_id, {}), f)

            # Hydrophobic residues
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_hydrophobic.pkl"), 'wb') as f:
                pickle.dump(data['hydrophobic_residues'].get(chain_id, {}), f)

            # Edge pairs
            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_edge_pairs.pkl"), 'wb') as f:
                pickle.dump(data['edge_indices'].get(chain_id, np.zeros((2, 0), dtype=int)), f)

            # Index file
            index = {
                'atoms': f"{pdb_id}_{chain_id}_atoms.pkl",
                'backbone': f"{pdb_id}_{chain_id}_backbone.pkl",
                'secondary_structure': f"{pdb_id}_{chain_id}_ss.pkl",
                'bond_lengths': f"{pdb_id}_{chain_id}_bond_lengths.pkl",
                'angles': f"{pdb_id}_{chain_id}_angles.pkl",
                'torsions': f"{pdb_id}_{chain_id}_torsions.pkl",
                'charges': f"{pdb_id}_{chain_id}_charges.pkl",
                'hydrophobic': f"{pdb_id}_{chain_id}_hydrophobic.pkl",
                'edge_pairs': f"{pdb_id}_{chain_id}_edge_pairs.pkl"
            }

            with open(os.path.join(chain_dir, f"{pdb_id}_{chain_id}_index.json"), 'w') as f:
                json.dump(index, f, indent=2)

def generate_synthetic_dataset(num_graphs=1000, nodes_per_graph=100, output_dir="synthetic_data"):
    """
    Generate a synthetic dataset of protein graphs and associated data.

    Parameters:
    -----------
    num_graphs : int
        Number of graphs to generate
    nodes_per_graph : int
        Number of nodes per graph
    output_dir : str
        Output directory
    """
    start_time = time.time()
    print(f"Generating {num_graphs} protein graphs with {nodes_per_graph} nodes each...")

    for i in tqdm(range(num_graphs)):
        # Use consistent seed for reproducibility
        seed = i
        pdb_id = f"synth{i:04d}"

        # Generate graph and data
        graph = generate_synthetic_protein_graph(nodes_per_graph, seed)
        data = generate_associated_data(graph, seed)

        # Save to disk
        save_graph_and_data(graph, data, output_dir, pdb_id)

    end_time = time.time()
    print(f"Dataset generation complete in {end_time - start_time:.1f} seconds")
    print(f"Data saved to {output_dir}")

if __name__ == "__main__":
    generate_synthetic_dataset(num_graphs=1000, nodes_per_graph=100)