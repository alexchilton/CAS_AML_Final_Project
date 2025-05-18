import os
import pickle
import pandas as pd
import numpy as np
import torch
from torch_geometric.data import Data, Dataset
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool, GATConv

# Load BRENDA data
path='C:/Users/laran/Documents/code/dataset/prot_selection/enzyme_structures_ph/enzyme_dataset_complete_ec_with_brenda_ph.csv'
brenda_data = pd.read_csv(path)
# Create a dictionary mapping PDB IDs to optimal pH values
pdb_to_ph = dict(zip(brenda_data['pdb_id'], brenda_data['ph_optimum']))

class ProteinGraphDataset(Dataset):
    def __init__(self, root, processed_dir, transform=None, pre_transform=None):
        # Initialize the parent class first
        super(ProteinGraphDataset, self).__init__(root, transform, pre_transform)
        
        # Store your directory in a private attribute
        self._custom_dir = processed_dir
        self.pdb_ids = []
        self.pdb_to_idx = {}
        
        print(f"Loading PDB data from: {processed_dir}")
        print(f"Number of PDB IDs in pH dictionary: {len(pdb_to_ph)}")
        
        # Find all processed PDB directories
        for pdb_dir in os.listdir(processed_dir):
            pdb_path = os.path.join(processed_dir, pdb_dir)
            if os.path.isdir(pdb_path):
                # Check if this PDB ID is in our pH dictionary
                if pdb_dir in pdb_to_ph:
                    # Check if the data file exists
                    data_file = os.path.join(pdb_path, f"{pdb_dir}_data.pkl")
                    graph_file = os.path.join(pdb_path, f"{pdb_dir}_graph.pkl")
                    
                    if os.path.exists(data_file) and os.path.exists(graph_file):
                        self.pdb_ids.append(pdb_dir)
                        self.pdb_to_idx[pdb_dir] = len(self.pdb_ids) - 1
                        #print(f"Added {pdb_dir} to dataset")
                    else:
                        # Check if we have chain-specific files
                        chain_files = [f for f in os.listdir(pdb_path) if f.startswith(f"{pdb_dir}_") and f.endswith("_data.pkl")]
                        if chain_files:
                            self.pdb_ids.append(pdb_dir)
                            self.pdb_to_idx[pdb_dir] = len(self.pdb_ids) - 1
                            #print(f"Added {pdb_dir} with chain files to dataset")
        
        print(f"Final dataset size: {len(self.pdb_ids)}")

    
    # Override the processed_dir property
    @property
    def processed_dir(self):
        return self._custom_dir

    @property
    def raw_file_names(self):
        return self.pdb_ids

    @property
    def processed_file_names(self):
        return [f'{pdb_id}.pt' for pdb_id in self.pdb_ids]

    def len(self):
        return len(self.pdb_ids)

    def get(self, idx):
        # Load from processed files if they exist
        processed_file = os.path.join(self.processed_dir, f'{self.pdb_ids[idx]}.pt')
        if os.path.exists(processed_file):
            return torch.load(processed_file)
        
        # Otherwise, create and process the data
        pdb_id = self.pdb_ids[idx]
        pdb_dir = os.path.join(self.processed_dir, pdb_id)
        
        # Check for main graph file
        graph_path = os.path.join(pdb_dir, f"{pdb_id}_graph.pkl")
        data_path = os.path.join(pdb_dir, f"{pdb_id}_data.pkl")
        
        if os.path.exists(graph_path) and os.path.exists(data_path):
            # Load main graph and data
            with open(graph_path, 'rb') as f:
                nx_graph = pickle.load(f)
            
            with open(data_path, 'rb') as f:
                pdb_data = pickle.load(f)
                
            return self.convert_to_pytorch_geometric(nx_graph, pdb_data, pdb_id)
        else:
            # Look for chain-specific files
            chain_ids = []
            for f in os.listdir(pdb_dir):
                # Look for files like "1A2V_A_index.json" to identify chains
                if f.startswith(f"{pdb_id}_") and f.endswith("_index.json"):
                    # Extract the chain ID (which is the character after the PDB ID)
                    parts = f.split('_')
                    if len(parts) >= 3:
                        chain_id = parts[1]  # Get the chain ID (like "A", "B", etc.)
                        chain_ids.append(chain_id)
            
            # Process each chain
            if chain_ids:
                # For simplicity, we'll use the first chain
                # In a real application, you might want to combine data from all chains
                chain_id = chain_ids[0]
                print(f"Using chain {chain_id} for {pdb_id}")
                
                # Construct paths for chain-specific files
                chain_graph_path = None  # We might not have separate graph files per chain
                
                # Collect chain-specific data
                chain_data = {}
                
                # Collect all relevant chain files
                for file_type in ["atoms", "backbone", "bond_lengths", "charges", 
                                "edge_pairs", "hydrophobic", "ss", "angles", "torsions"]:
                    file_path = os.path.join(pdb_dir, f"{pdb_id}_{chain_id}_{file_type}.pkl")
                    if os.path.exists(file_path):
                        with open(file_path, 'rb') as f:
                            chain_data[file_type] = pickle.load(f)
                
                # We need to check if there's a main graph file we can use
                graph_path = os.path.join(pdb_dir, f"{pdb_id}_graph.pkl")
                if os.path.exists(graph_path):
                    with open(graph_path, 'rb') as f:
                        nx_graph = pickle.load(f)
                    
                    # Convert chain data to appropriate format
                    pdb_data = {
                        'pdb_id': pdb_id,
                        'chain_id': chain_id,
                        'residues_by_chain': {chain_id: chain_data.get('atoms', [])},
                        'secondary_structure': {chain_id: chain_data.get('ss', {})},
                        'backbone_atoms': {chain_id: chain_data.get('backbone', {})},
                        'bond_lengths': {chain_id: chain_data.get('bond_lengths', {})},
                        'charges': {chain_id: chain_data.get('charges', {})}
                    }
                    
                    return self.convert_to_pytorch_geometric(nx_graph, pdb_data, pdb_id)
                else:
                    # If we don't have a graph file, we need to construct the graph
                    # This is more complex and depends on your preprocessing
                    # For now, raise an error
                    raise FileNotFoundError(f"Could not find graph file for {pdb_id}")
            
            # If we reach here, we couldn't process the data
            raise FileNotFoundError(f"Could not find data for {pdb_id}")

    def convert_to_pytorch_geometric(self, nx_graph, pdb_data, pdb_id):
        # Extract node features
        x = []
        node_mapping = {}  # Map original node IDs to consecutive indices
        
        for i, (node, attrs) in enumerate(nx_graph.nodes(data=True)):
            node_mapping[node] = i
            
            # Create node feature vector
            features = [
                # One-hot encode residue type
                self.one_hot_residue(attrs.get('residue_name', 'UNK')),
                
                # Secondary structure (one-hot)
                self.one_hot_ss(attrs.get('ss', '?')),
                
                # Add hydrophobicity (binary)
                1.0 if attrs.get('residue_name') in ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'] else 0.0,
                
                # Add charge-related features
                self.get_residue_charge(attrs.get('residue_name', 'UNK')),
            ]
            
            # Flatten feature vector
            x.append(np.concatenate(features))
        
        # Create edge index and edge features
        edge_index = []
        edge_attr = []
        
        for u, v, edge_data in nx_graph.edges(data=True):
            # Skip if nodes not in mapping (should not happen)
            if u not in node_mapping or v not in node_mapping:
                continue
                
            # Add edge in both directions (undirected graph)
            edge_index.append([node_mapping[u], node_mapping[v]])
            edge_index.append([node_mapping[v], node_mapping[u]])
            
            # Edge features
            edge_features = [
                1.0 if 'peptide_bond' in edge_data.get('kind', {}) else 0.0,
                1.0 if 'contact' in edge_data.get('kind', {}) else 0.0,
                edge_data.get('distance', 0.0) if edge_data.get('distance') is not None else 0.0,
            ]
            
            # Add edge features in both directions
            edge_attr.append(edge_features)
            edge_attr.append(edge_features)
        
        # Convert to PyTorch tensors
        x = torch.tensor(np.array(x), dtype=torch.float)
        edge_index = torch.tensor(np.array(edge_index).T, dtype=torch.long)
        edge_attr = torch.tensor(np.array(edge_attr), dtype=torch.float)
        
        # Get the target pH value
        y = torch.tensor([pdb_to_ph.get(pdb_id, 7.0)], dtype=torch.float)
        
        # Create PyTorch Geometric Data object
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y)
        data.pdb_id = pdb_id
        
        # Save processed data
        torch.save(data, os.path.join(self.processed_dir, f'{pdb_id}.pt'))
        
        return data
    
    def one_hot_residue(self, residue_name):
        # 20 standard amino acids + unknown
        amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
                      'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'UNK']
        
        if residue_name not in amino_acids:
            residue_name = 'UNK'
            
        one_hot = [0] * len(amino_acids)
        one_hot[amino_acids.index(residue_name)] = 1
        
        return one_hot
    
    def one_hot_ss(self, ss_type):
        # Secondary structure types: H (alpha helix), E (beta strand), C (coil), etc.
        ss_types = ['H', 'E', 'C', 'B', 'G', 'I', 'T', 'S', '?']
        
        if ss_type not in ss_types:
            ss_type = '?'
            
        one_hot = [0] * len(ss_types)
        one_hot[ss_types.index(ss_type)] = 1
        
        return one_hot
    
    def get_residue_charge(self, residue_name):
        # Assign charge values to amino acids based on their properties at physiological pH
        charge_map = {
            'ARG': 1.0,  # Positive
            'LYS': 1.0,  # Positive
            'HIS': 0.1,  # Slightly positive at physiological pH
            'ASP': -1.0, # Negative
            'GLU': -1.0, # Negative
            'CYS': -0.1, # Slightly negative when deprotonated
            'TYR': -0.1, # Slightly negative when deprotonated
        }
        
        return [charge_map.get(residue_name, 0.0)]