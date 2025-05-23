import os
import pickle
import json
import numpy as np
import networkx as nx
import torch
from torch_geometric.data import Data
import traceback
from tqdm import tqdm
import gc

class NanobodyGraphBuilder:
    """
    A class for integrating physical properties into NetworkX graphs
    and converting them to PyTorch Geometric format for PINN GNN training.
    """
    
    def __init__(self, data_dir):
        """
        Initialize the graph builder.
        
        Parameters:
        -----------
        data_dir : str
            Directory containing processed nanobody data
        """
        self.data_dir = data_dir
        self.processed_dir = os.path.join(data_dir, "processed_data")
        
        # Feature normalization parameters
        self.feature_stats = {
            'bond_lengths': {'mean': None, 'std': None},
            'angles': {'mean': None, 'std': None},
            'torsions': {'mean': None, 'std': None},
            'charges': {'mean': None, 'std': None}
        }
        
        # Secondary structure encoding
        self.ss_encoding = {
            'H': [1, 0, 0, 0, 0, 0, 0, 0],  # Alpha helix
            'B': [0, 1, 0, 0, 0, 0, 0, 0],  # Beta bridge
            'E': [0, 0, 1, 0, 0, 0, 0, 0],  # Extended strand
            'G': [0, 0, 0, 1, 0, 0, 0, 0],  # 3-10 helix
            'I': [0, 0, 0, 0, 1, 0, 0, 0],  # Pi helix
            'T': [0, 0, 0, 0, 0, 1, 0, 0],  # Turn
            'S': [0, 0, 0, 0, 0, 0, 1, 0],  # Bend
            '?': [0, 0, 0, 0, 0, 0, 0, 1],  # Unknown
            ' ': [0, 0, 0, 0, 0, 0, 0, 1]   # Unknown
        }
        
        # Amino acid encoding (one-hot)
        self.aa_list = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK'
        ]
        self.aa_encoding = {aa: i for i, aa in enumerate(self.aa_list)}
    
    def get_processed_nanobodies(self):
        """
        Get a list of all processed nanobody directories.
        
        Returns:
        --------
        list
            List of nanobody IDs
        """
        if not os.path.exists(self.processed_dir):
            raise ValueError(f"Processed data directory not found: {self.processed_dir}")
        
        nanobody_dirs = [d for d in os.listdir(self.processed_dir) 
                        if os.path.isdir(os.path.join(self.processed_dir, d))]
        
        print(f"Found {len(nanobody_dirs)} processed nanobodies")
        return nanobody_dirs
    
    def load_nanobody_data(self, nanobody_id):
        """
        Load all data for a single nanobody, including chain-specific files.
        
        Parameters:
        -----------
        nanobody_id : str
            ID of the nanobody to load
            
        Returns:
        --------
        tuple
            (graph, data_dict) containing the NetworkX graph and data dictionary
        """
        nanobody_dir = os.path.join(self.processed_dir, nanobody_id)
        
        # Load main graph
        graph_path = os.path.join(nanobody_dir, f"{nanobody_id}_graph.pkl")
        if not os.path.exists(graph_path):
            raise ValueError(f"Graph file not found: {graph_path}")
        
        with open(graph_path, 'rb') as f:
            graph = pickle.load(f)
        
        # Load main data dictionary
        data_path = os.path.join(nanobody_dir, f"{nanobody_id}_data.pkl")
        if not os.path.exists(data_path):
            raise ValueError(f"Data file not found: {data_path}")
        
        with open(data_path, 'rb') as f:
            data_dict = pickle.load(f)
        
        # Enhance data_dict with chain-specific data
        chains = list(data_dict.get('residues_by_chain', {}).keys())
        
        for chain_id in chains:
            # Load chain index file to find all available chain data
            index_path = os.path.join(nanobody_dir, f"{nanobody_id}_{chain_id}_index.json")
            
            if os.path.exists(index_path):
                # Load the chain index
                with open(index_path, 'r') as f:
                    chain_index = json.load(f)
                
                # Initialize chain data in the main dictionary if needed
                for key in ['backbone_atoms', 'secondary_structure', 'bond_lengths', 
                        'angles', 'torsions', 'charges', 'hydrophobic_residues', 'edge_indices']:
                    if key not in data_dict:
                        data_dict[key] = {}
                
                # Load and add chain-specific data
                if 'backbone' in chain_index and os.path.exists(chain_index['backbone']):
                    with open(chain_index['backbone'], 'rb') as f:
                        data_dict['backbone_atoms'][chain_id] = pickle.load(f)
                
                if 'secondary_structure' in chain_index and os.path.exists(chain_index['secondary_structure']):
                    with open(chain_index['secondary_structure'], 'rb') as f:
                        data_dict['secondary_structure'][chain_id] = pickle.load(f)
                
                if 'bond_lengths' in chain_index and os.path.exists(chain_index['bond_lengths']):
                    with open(chain_index['bond_lengths'], 'rb') as f:
                        data_dict['bond_lengths'][chain_id] = pickle.load(f)
                
                if 'angles' in chain_index and os.path.exists(chain_index['angles']):
                    with open(chain_index['angles'], 'rb') as f:
                        data_dict['angles'][chain_id] = pickle.load(f)
                
                if 'torsions' in chain_index and os.path.exists(chain_index['torsions']):
                    with open(chain_index['torsions'], 'rb') as f:
                        data_dict['torsions'][chain_id] = pickle.load(f)
                
                if 'charges' in chain_index and os.path.exists(chain_index['charges']):
                    with open(chain_index['charges'], 'rb') as f:
                        data_dict['charges'][chain_id] = pickle.load(f)
                
                if 'hydrophobic' in chain_index and os.path.exists(chain_index['hydrophobic']):
                    with open(chain_index['hydrophobic'], 'rb') as f:
                        data_dict['hydrophobic_residues'][chain_id] = pickle.load(f)
                
                if 'edge_pairs' in chain_index and os.path.exists(chain_index['edge_pairs']):
                    with open(chain_index['edge_pairs'], 'rb') as f:
                        data_dict['edge_indices'][chain_id] = pickle.load(f)
                
                # Optionally load atoms if needed
                if 'atoms' in chain_index and os.path.exists(chain_index['atoms']):
                    if 'atoms' not in data_dict:
                        data_dict['atoms'] = {}
                    with open(chain_index['atoms'], 'rb') as f:
                        data_dict['atoms'][chain_id] = pickle.load(f)
        
        return graph, data_dict


    def integrate_properties_into_graph(self, graph, data_dict):
        """
        Integrate all physical properties from the data dictionary into the graph.
        Specifically designed to handle the nanobody data format with chain-specific files.
        
        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to update
        data_dict : dict
            Dictionary of physical properties
            
        Returns:
        --------
        networkx.Graph
            Updated graph with integrated properties
        """
        # Make a copy of the graph to avoid modifying the original
        enhanced_graph = graph.copy()
        
        # Extract chain IDs
        chain_ids = list(data_dict.get('residues_by_chain', {}).keys())
        
        # Add residue-level properties to nodes
        for node, data in enhanced_graph.nodes(data=True):
            # Extract chain and residue information
            chain_id, res_id, res_name = None, None, None
            
            if ':' in node:  # Format like 'H:VAL:1'
                parts = node.split(':')
                if len(parts) >= 3:
                    chain_id = parts[0]
                    res_name = parts[1]
                    try:
                        res_id = int(parts[2])
                    except ValueError:
                        res_id = parts[2]
            else:
                # Try to extract from attributes
                chain_id = data.get('chain_id', chain_ids[0] if chain_ids else '')
                res_id = data.get('residue_number')
                res_name = data.get('residue_name', 'UNK')
            
            # Skip if missing critical information
            if chain_id is None or res_id is None:
                continue
            
            # 1. Add amino acid type as one-hot encoding
            aa_index = self.aa_encoding.get(res_name, self.aa_encoding.get('UNK', len(self.aa_list) - 1))
            enhanced_graph.nodes[node]['aa_type'] = aa_index
            
            # 2. Add secondary structure (if available)
            if 'secondary_structure' in data_dict and chain_id in data_dict['secondary_structure']:
                ss = data_dict['secondary_structure'][chain_id].get(res_id, '?')
                if isinstance(ss, dict) and 'ss' in ss:  # Handle different formats
                    ss = ss['ss']
                enhanced_graph.nodes[node]['ss'] = ss
                
                # Add as one-hot encoding
                if ss in self.ss_encoding:
                    enhanced_graph.nodes[node]['ss_encoding'] = self.ss_encoding[ss]
                else:
                    enhanced_graph.nodes[node]['ss_encoding'] = self.ss_encoding['?']
            else:
                # Default secondary structure
                enhanced_graph.nodes[node]['ss'] = '?'
                enhanced_graph.nodes[node]['ss_encoding'] = self.ss_encoding['?']
            
            # 3. Add hydrophobicity information
            if 'hydrophobic_residues' in data_dict and chain_id in data_dict['hydrophobic_residues']:
                hydro_set = data_dict['hydrophobic_residues'][chain_id]
                # Handle both set and list formats
                if isinstance(hydro_set, (set, list)):
                    is_hydrophobic = res_id in hydro_set
                else:
                    # Try direct dictionary access
                    is_hydrophobic = hydro_set.get(res_id, False)
                enhanced_graph.nodes[node]['is_hydrophobic'] = float(is_hydrophobic)
            else:
                enhanced_graph.nodes[node]['is_hydrophobic'] = 0.0  # Default: not hydrophobic
            
            # 4. Add backbone atom positions (if available)
            if 'backbone_atoms' in data_dict and chain_id in data_dict['backbone_atoms']:
                if res_id in data_dict['backbone_atoms'][chain_id]:
                    backbone = data_dict['backbone_atoms'][chain_id][res_id]
                    for atom_name, position in backbone.items():
                        enhanced_graph.nodes[node][f'backbone_{atom_name}'] = position
                    
                    # Mark that this node has backbone atoms
                    enhanced_graph.nodes[node]['has_backbone'] = True
                    
                    # Set coordinates (for visualization or spatial calculations)
                    if 'CA' in backbone:
                        enhanced_graph.nodes[node]['coords'] = backbone['CA']
                else:
                    enhanced_graph.nodes[node]['has_backbone'] = False
            else:
                enhanced_graph.nodes[node]['has_backbone'] = False
            
            # 5. Add charge information (if available)
            if 'charges' in data_dict and chain_id in data_dict['charges']:
                node_charges = {}
                
                # If charges are per atom
                atom_charges = data_dict['charges'][chain_id]
                if 'atoms' in data_dict and chain_id in data_dict['atoms']:
                    atoms = data_dict['atoms'][chain_id]
                    
                    # Find atoms for this residue
                    for i, atom in enumerate(atoms):
                        atom_chain, atom_res, atom_name = atom[0], atom[1], atom[3]  # Unpack based on your atom format
                        if atom_chain == chain_id and atom_res == res_id:
                            # Check if this atom index is in charges
                            if i in atom_charges:
                                node_charges[atom_name] = atom_charges[i]
                
                # If node_charges is empty, try direct access by residue ID
                if not node_charges and res_id in atom_charges:
                    if isinstance(atom_charges[res_id], dict):
                        node_charges = atom_charges[res_id]
                    else:
                        node_charges = {'mean': atom_charges[res_id]}
                
                # Set mean charge
                if node_charges:
                    enhanced_graph.nodes[node]['mean_charge'] = sum(node_charges.values()) / len(node_charges)
                else:
                    enhanced_graph.nodes[node]['mean_charge'] = 0.0
            else:
                enhanced_graph.nodes[node]['mean_charge'] = 0.0
            
            # 6. Add angle and torsion statistics (if available)
            if 'angles' in data_dict and chain_id in data_dict['angles']:
                angles_dict = data_dict['angles'][chain_id]
                # Find angles where this residue is the central atom
                res_angles = []
                
                for angle_key, angle_val in angles_dict.items():
                    # Check if this is a tuple format like (i, j, k)
                    if isinstance(angle_key, tuple) and len(angle_key) == 3:
                        # If this residue is the middle of the angle
                        if angle_key[1] == res_id:
                            res_angles.append(angle_val)
                
                if res_angles:
                    enhanced_graph.nodes[node]['mean_angle'] = sum(res_angles) / len(res_angles)
                    enhanced_graph.nodes[node]['min_angle'] = min(res_angles)
                    enhanced_graph.nodes[node]['max_angle'] = max(res_angles)
            
            if 'torsions' in data_dict and chain_id in data_dict['torsions']:
                torsions_dict = data_dict['torsions'][chain_id]
                # Find torsions where this residue is part of the sequence
                res_torsions = []
                
                for torsion_key, torsion_val in torsions_dict.items():
                    # Check if this is a tuple format like (i, j, k, l)
                    if isinstance(torsion_key, tuple) and len(torsion_key) == 4:
                        # If this residue is in the torsion
                        if res_id in torsion_key:
                            res_torsions.append(torsion_val)
                
                if res_torsions:
                    enhanced_graph.nodes[node]['mean_torsion'] = sum(res_torsions) / len(res_torsions)
                    enhanced_graph.nodes[node]['min_torsion'] = min(res_torsions)
                    enhanced_graph.nodes[node]['max_torsion'] = max(res_torsions)
        
        # Add edge properties
        for edge in enhanced_graph.edges():
            u, v = edge
            u_chain, u_res = None, None
            v_chain, v_res = None, None
            
            # Extract residue info from nodes
            if ':' in u:
                parts = u.split(':')
                if len(parts) >= 3:
                    u_chain = parts[0]
                    try:
                        u_res = int(parts[2])
                    except ValueError:
                        u_res = parts[2]
            
            if ':' in v:
                parts = v.split(':')
                if len(parts) >= 3:
                    v_chain = parts[0]
                    try:
                        v_res = int(parts[2])
                    except ValueError:
                        v_res = parts[2]
            
            # Skip if we couldn't extract needed information
            if u_chain is None or v_chain is None or u_res is None or v_res is None:
                continue
            
            # 1. Add bond length information (if available)
            if 'bond_lengths' in data_dict:
                # Only consider if both residues are in the same chain
                if u_chain == v_chain and u_chain in data_dict['bond_lengths']:
                    bond_lengths = data_dict['bond_lengths'][u_chain]
                    
                    # Find bond lengths between atoms of these residues
                    for (i, j), length in bond_lengths.items():
                        # We need to determine if atoms i and j belong to residues u_res and v_res
                        # This requires atom information
                        if 'atoms' in data_dict and u_chain in data_dict['atoms']:
                            atoms = data_dict['atoms'][u_chain]
                            
                            # Only proceed if the atom indices are in range
                            if i < len(atoms) and j < len(atoms):
                                i_chain, i_res = atoms[i][0], atoms[i][1]
                                j_chain, j_res = atoms[j][0], atoms[j][1]
                                
                                # Check if these atoms connect our residues
                                if ((i_res == u_res and j_res == v_res) or
                                    (i_res == v_res and j_res == u_res)):
                                    enhanced_graph.edges[u, v]['bond_length'] = length
                                    break
            
            # 2. Tag sequential peptide bonds
            if u_chain == v_chain:
                is_sequential = abs(u_res - v_res) == 1
                enhanced_graph.edges[u, v]['is_peptide_bond'] = float(is_sequential)
            else:
                enhanced_graph.edges[u, v]['is_peptide_bond'] = 0.0
            
            # 3. Flag inter-chain edges
            enhanced_graph.edges[u, v]['is_interchain'] = float(u_chain != v_chain)
        
        return enhanced_graph

    
    


    
    def _extract_node_info(self, node, node_data):
        """Extract chain ID and residue ID from a node."""
        if ':' in node:  # Format like 'H:VAL:1'
            parts = node.split(':')
            chain_id = parts[0]
            try:
                res_id = int(parts[2])
                return chain_id, res_id
            except (ValueError, IndexError):
                return None
        else:
            # Try to extract from attributes
            chain_id = node_data.get('chain_id')
            res_id = node_data.get('residue_number')
            if chain_id is not None and res_id is not None:
                return chain_id, res_id
        
        return None
    
    def convert_to_pytorch_geometric(self, graph, nanobody_id, pH=None, ionic_strength=None):
        """
        Convert a NetworkX graph to PyTorch Geometric format.
        
        Parameters:
        -----------
        graph : networkx.Graph
            Enhanced NetworkX graph with integrated properties
        nanobody_id : str
            ID of the nanobody
        pH : float, optional
            pH value (for labeled data)
        ionic_strength : float, optional
            Ionic strength value (for labeled data)
            
        Returns:
        --------
        torch_geometric.data.Data
            PyTorch Geometric graph
        """
        # 1. Create node feature matrix
        node_list = list(graph.nodes())
        num_nodes = len(node_list)
        
        # Define features to extract
        node_features = []
        
        for node in node_list:
            node_data = graph.nodes[node]
            
            # Create feature vector for this node
            features = []
            
            # 1. One-hot encoded amino acid type
            aa_type = node_data.get('aa_type', self.aa_encoding['UNK'])
            aa_one_hot = [0] * len(self.aa_list)
            aa_one_hot[aa_type] = 1
            features.extend(aa_one_hot)
            
            # 2. Secondary structure encoding
            ss_encoding = node_data.get('ss_encoding', self.ss_encoding['?'])
            features.extend(ss_encoding)
            
            # 3. Hydrophobicity
            is_hydrophobic = node_data.get('is_hydrophobic', 0.0)
            features.append(is_hydrophobic)
            
            # 4. Charge information
            mean_charge = node_data.get('mean_charge', 0.0)
            features.append(mean_charge)
            
            # 5. Backbone atom positions (if used as features rather than positions)
            # Note: Usually, you would use these as the node positions, not features
            # Here we're just adding placeholders
            has_backbone = node_data.get('has_backbone', False)
            features.append(float(has_backbone))
            
            # Add to list of node features
            node_features.append(features)
        
        # 2. Create edge index
        edge_index = []
        edge_attr = []
        
        for u, v, data in graph.edges(data=True):
            # Get node indices
            u_idx = node_list.index(u)
            v_idx = node_list.index(v)
            
            # Add edges in both directions (undirected graph)
            edge_index.append([u_idx, v_idx])
            edge_index.append([v_idx, u_idx])
            
            # Create edge features
            edge_features = []
            
            # 1. Bond length
            bond_length = data.get('bond_length', 0.0)
            edge_features.append(bond_length)
            
            # 2. Is peptide bond
            is_peptide = data.get('is_peptide_bond', 0.0)
            edge_features.append(is_peptide)
            
            # 3. Has angle info
            has_angle = data.get('has_angle_info', False)
            edge_features.append(float(has_angle))
            
            # Add same features for reverse edge
            edge_attr.append(edge_features)
            edge_attr.append(edge_features)
        
        # 3. Create node positions
        pos = []
        for node in node_list:
            node_data = graph.nodes[node]
            
            # Try to get CA atom position
            if 'backbone_CA' in node_data:
                position = node_data['backbone_CA']
            elif 'CA' in node_data:
                position = node_data['CA']
            elif 'coords' in node_data:
                position = node_data['coords']
            else:
                # Fallback position
                position = np.zeros(3)
            
            pos.append(position)
        
        # 4. Create PyTorch tensors
        x = torch.tensor(node_features, dtype=torch.float)
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)
        pos = torch.tensor(pos, dtype=torch.float)
        
        # 5. Add labels if provided
        y = None
        if pH is not None and ionic_strength is not None:
            y = torch.tensor([pH, ionic_strength], dtype=torch.float)
        
        # 6. Create PyTorch Geometric data object
        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
            pos=pos,
            y=y,
            nanobody_id=nanobody_id
        )
        
        return data
    
    def process_all_nanobodies(self, output_dir=None, include_labels=False, label_file=None):
        """
        Process all nanobodies and convert them to PyTorch Geometric format.
        
        Parameters:
        -----------
        output_dir : str, optional
            Directory to save processed graphs
        include_labels : bool
            Whether to include pH and ionic strength labels
        label_file : str, optional
            Path to JSON file containing labels for nanobodies
            
        Returns:
        --------
        list
            List of PyTorch Geometric Data objects
        """
        if output_dir is None:
            output_dir = os.path.join(self.data_dir, "pytorch_geometric")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Load labels if provided
        labels = {}
        if include_labels and label_file and os.path.exists(label_file):
            with open(label_file, 'r') as f:
                labels = json.load(f)
        
        # Get list of nanobodies
        nanobody_ids = self.get_processed_nanobodies()
        
        # Process each nanobody
        geometric_graphs = []
        
        for nanobody_id in tqdm(nanobody_ids, desc="Processing nanobodies"):
            try:
                # Load data
                graph, data_dict = self.load_nanobody_data(nanobody_id)
                
                # Integrate properties
                enhanced_graph = self.integrate_properties_into_graph(graph, data_dict)
                
                # Get labels if available
                pH = None
                ionic_strength = None
                
                if include_labels and nanobody_id in labels:
                    pH = labels[nanobody_id].get('optimal_pH')
                    ionic_strength = labels[nanobody_id].get('optimal_ionic_strength')
                
                # Convert to PyTorch Geometric
                geo_data = self.convert_to_pytorch_geometric(
                    enhanced_graph, nanobody_id, pH, ionic_strength
                )
                
                # Save to file
                output_file = os.path.join(output_dir, f"{nanobody_id}.pt")
                torch.save(geo_data, output_file)
                
                # Add to list
                geometric_graphs.append(geo_data)
                
                # Clean up to avoid memory issues
                del graph, data_dict, enhanced_graph, geo_data
                gc.collect()
                
            except Exception as e:
                print(f"Error processing {nanobody_id}: {str(e)}")
                traceback.print_exc()
        
        return geometric_graphs
    
    def create_dataset(self, train_split=0.8, val_split=0.1, test_split=0.1, seed=42):
        """
        Create train/val/test splits of the dataset.
        
        Parameters:
        -----------
        train_split : float
            Proportion of data for training
        val_split : float
            Proportion of data for validation
        test_split : float
            Proportion of data for testing
        seed : int
            Random seed for reproducibility
            
        Returns:
        --------
        dict
            Dictionary with 'train', 'val', and 'test' keys
        """
        # Set random seed
        np.random.seed(seed)
        
        # Get list of all processed files
        dataset_dir = os.path.join(self.data_dir, "pytorch_geometric")
        all_files = [f for f in os.listdir(dataset_dir) if f.endswith('.pt')]
        
        # Shuffle files
        np.random.shuffle(all_files)
        
        # Calculate split indices
        n_train = int(len(all_files) * train_split)
        n_val = int(len(all_files) * val_split)
        
        # Split files
        train_files = all_files[:n_train]
        val_files = all_files[n_train:n_train+n_val]
        test_files = all_files[n_train+n_val:]
        
        # Create dataset splits
        splits = {
            'train': [os.path.join(dataset_dir, f) for f in train_files],
            'val': [os.path.join(dataset_dir, f) for f in val_files],
            'test': [os.path.join(dataset_dir, f) for f in test_files]
        }
        
        # Save split information
        split_info = {
            'train': [os.path.splitext(f)[0] for f in train_files],
            'val': [os.path.splitext(f)[0] for f in val_files],
            'test': [os.path.splitext(f)[0] for f in test_files]
        }
        
        with open(os.path.join(dataset_dir, 'splits.json'), 'w') as f:
            json.dump(split_info, f, indent=2)
        
        return splits

# Example usage
if __name__ == "__main__":
    # Set the path to your nanobody data directory
    data_dir = "nanos_networkx"
    
    # Create graph builder
    builder = NanobodyGraphBuilder(data_dir)
    
    # Process all nanobodies
    geometric_graphs = builder.process_all_nanobodies(
        include_labels=True,
        label_file="nanobody_optimal_conditions.json"  # Optional: path to file with pH/ionic strength labels
    )
    
    # Create dataset splits
    splits = builder.create_dataset(train_split=0.8, val_split=0.1, test_split=0.1)
    
    print(f"Created dataset with {len(splits['train'])} training, {len(splits['val'])} validation, and {len(splits['test'])} test samples")