import os
import pickle
import json
import numpy as np
import torch
from torch_geometric.data import Data, InMemoryDataset
import networkx as nx
from tqdm import tqdm
import time


class ProteinGraphDataset(object):
    """
    Custom dataset for loading synthetic protein graphs
    (Not inheriting from PyTorch Geometric Dataset to avoid conflicts)
    """
    def __init__(self, root_dir, max_num_graphs=None):
        """
        Args:
            root_dir: Directory containing the synthetic protein data
            max_num_graphs: Maximum number of graphs to load (None for all)
        """
        self.root_dir = root_dir
        self.max_num_graphs = max_num_graphs

        # Find all PDB IDs
        self.pdb_ids = []
        for entry in os.listdir(root_dir):
            if os.path.isdir(os.path.join(root_dir, entry)) and entry.startswith("synth"):
                self.pdb_ids.append(entry)

        # Limit to specified number if needed
        if max_num_graphs is not None:
            self.pdb_ids = self.pdb_ids[:max_num_graphs]

        print(f"Found {len(self.pdb_ids)} protein graphs in {root_dir}")

    def __len__(self):
        return len(self.pdb_ids)

    def __getitem__(self, idx):
        """Load a single graph by index"""
        pdb_id = self.pdb_ids[idx]
        return self._load_graph(pdb_id)

    def _load_graph(self, pdb_id):
        """Load a graph and convert to PyTorch Geometric Data"""
        pdb_dir = os.path.join(self.root_dir, pdb_id)

        # Load NetworkX graph
        graph_path = os.path.join(pdb_dir, f"{pdb_id}_graph.pkl")
        with open(graph_path, 'rb') as f:
            nx_graph = pickle.load(f)

        # Load data dictionary
        data_path = os.path.join(pdb_dir, f"{pdb_id}_data.pkl")
        with open(data_path, 'rb') as f:
            data_dict = pickle.load(f)

        # Convert to PyTorch Geometric
        return self._convert_to_pytorch_geometric(nx_graph, data_dict, pdb_id)

    def _convert_to_pytorch_geometric(self, nx_graph, data_dict, pdb_id):
        """Convert NetworkX graph to PyTorch Geometric Data"""
        # Create node mapping
        node_mapping = {node: i for i, node in enumerate(nx_graph.nodes())}

        # Get edge indices
        edge_index = []
        edge_attr = []

        for u, v, data in nx_graph.edges(data=True):
            # Add edge in both directions (PyTorch Geometric convention)
            edge_index.append([node_mapping[u], node_mapping[v]])
            edge_index.append([node_mapping[v], node_mapping[u]])

            # Edge features (distance)
            distance = data.get('distance', 1.0)
            edge_attr.append([distance])
            edge_attr.append([distance])

        # Convert to tensors
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)

        # Create node features
        node_features = []

        for node in nx_graph.nodes():
            data = nx_graph.nodes[node]
            features = []

            # Basic features always included
            # 1. One-hot encode amino acid
            amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
            res_name = data.get('residue_name', 'UNK')
            one_hot = [1.0 if res_name == aa else 0.0 for aa in amino_acids]
            features.extend(one_hot)

            # 2. Secondary structure
            ss_types = ['H', 'E', 'B', 'G', 'T', 'S', '?']
            ss = data.get('ss', '?')
            ss_one_hot = [1.0 if ss == ss_type else 0.0 for ss_type in ss_types]
            features.extend(ss_one_hot)

            # 3. Backbone flag
            features.append(1.0 if data.get('has_backbone', False) else 0.0)

            # 4. 3D coordinates
            if 'coords' in data:
                coords = data['coords']
                # Normalize coordinates
                norm_coords = np.tanh(coords / 10.0)  # Scale to approximately [-1, 1]
                features.extend(norm_coords)
            else:
                features.extend([0.0, 0.0, 0.0])

            # 5. Additional features from data dictionary
            chain_id = data.get('chain_id', '')
            res_id = data.get('residue_number', 0)

            # Hydrophobicity
            is_hydrophobic = False
            if 'hydrophobic_residues' in data_dict and chain_id in data_dict['hydrophobic_residues']:
                is_hydrophobic = res_id in data_dict['hydrophobic_residues'][chain_id]
            features.append(1.0 if is_hydrophobic else 0.0)

            # Average bond length
            avg_bond_length = 0.0
            count = 0
            if 'bond_lengths' in data_dict and chain_id in data_dict['bond_lengths']:
                for (i, j), length in data_dict['bond_lengths'][chain_id].items():
                    count += 1
                    avg_bond_length += length

            if count > 0:
                avg_bond_length /= count
            features.append(avg_bond_length)

            node_features.append(features)

        # Convert to tensor
        x = torch.tensor(node_features, dtype=torch.float)

        # Create PyTorch Geometric Data object
        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
            num_nodes=len(node_features)
        )

        # Add metadata
        data.pdb_id = pdb_id

        return data


class ProteinGraphInMemoryDataset(InMemoryDataset):
    """
    PyTorch Geometric InMemoryDataset for protein graphs
    """
    def __init__(self, root_dir, transform=None, pre_transform=None, max_num_graphs=None):
        self.root_dir = root_dir
        self.max_num_graphs = max_num_graphs

        # Load data directly without using the standard PyTorch Geometric pattern
        # as we already have the data processed
        self.data_list = self._load_all_data()

        super(ProteinGraphInMemoryDataset, self).__init__(root_dir, transform, pre_transform)

    def _load_all_data(self):
        """Load all graphs into memory"""
        # Use the custom loader
        dataset = ProteinGraphDataset(self.root_dir, self.max_num_graphs)

        # Load all graphs
        print("Loading all graphs into memory...")
        all_data = []
        for i in tqdm(range(len(dataset))):
            all_data.append(dataset[i])

        return all_data

    def len(self):
        return len(self.data_list)

    def get(self, idx):
        return self.data_list[idx]


def load_all_graphs(data_dir, max_num=None, batch_size=32):
    """Load all graphs into memory and track memory usage"""
    # Print memory usage summary before loading
    try:
        import psutil
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / (1024 * 1024)  # MB
        print(f"Memory usage before loading: {mem_before:.2f} MB")
    except ImportError:
        print("psutil not available for memory tracking")

    # Load all graphs
    start_time = time.time()
    dataset = ProteinGraphDataset(data_dir, max_num_graphs=max_num)

    all_data = []
    for i in tqdm(range(len(dataset))):
        all_data.append(dataset[i])

    end_time = time.time()
    print(f"Loading time: {end_time - start_time:.2f} seconds")

    # Print memory usage after loading
    try:
        import psutil
        process = psutil.Process(os.getpid())
        mem_after = process.memory_info().rss / (1024 * 1024)  # MB
        print(f"Memory usage after loading: {mem_after:.2f} MB")
        print(f"Memory increase: {mem_after - mem_before:.2f} MB")
    except ImportError:
        pass

    # Summarize dataset
    total_nodes = sum(data.num_nodes for data in all_data)
    total_edges = sum(data.edge_index.size(1) for data in all_data)
    avg_nodes = total_nodes / len(all_data)
    avg_edges = total_edges / len(all_data)

    print(f"Loaded {len(all_data)} graphs with:")
    print(f"- Total nodes: {total_nodes}")
    print(f"- Total edges: {total_edges}")
    print(f"- Average nodes per graph: {avg_nodes:.2f}")
    print(f"- Average edges per graph: {avg_edges:.2f}")
    print(f"- Feature dimensions: {all_data[0].x.size(1)}")

    return dataset, all_data


def test_batch_loading(all_data, batch_size=32):
    """Test loading and batching graphs"""
    from torch_geometric.loader import DataLoader

    # Create in-memory dataset
    class SimpleDataset(torch.utils.data.Dataset):
        def __init__(self, data_list):
            self.data_list = data_list

        def __len__(self):
            return len(self.data_list)

        def __getitem__(self, idx):
            return self.data_list[idx]

    simple_dataset = SimpleDataset(all_data)
    loader = DataLoader(simple_dataset, batch_size=batch_size, shuffle=True)

    print(f"\nTesting batch loading with batch_size={batch_size}...")
    start_time = time.time()

    for i, batch in enumerate(loader):
        # Print info for first batch
        if i == 0:
            print(f"Batch info:")
            print(f"- Nodes: {batch.x.size(0)}")
            print(f"- Edges: {batch.edge_index.size(1)}")
            print(f"- Node features: {batch.x.size(1)}")
            print(f"- Number of graphs in batch: {batch.num_graphs}")

        # Process a few batches then break
        if i >= 3:
            break

    end_time = time.time()
    print(f"Batch loading test completed in {end_time - start_time:.3f} seconds")

    # Test a forward pass through a simple GNN
    try:
        from torch_geometric.nn import GCNConv, global_mean_pool
        import torch.nn.functional as F

        # Simple GNN model
        class SimpleGNN(torch.nn.Module):
            def __init__(self, in_channels):
                super(SimpleGNN, self).__init__()
                self.conv1 = GCNConv(in_channels, 64)
                self.conv2 = GCNConv(64, 64)
                self.fc = torch.nn.Linear(64, 2)  # Binary classification

            def forward(self, x, edge_index, batch):
                # GNN layers
                x = self.conv1(x, edge_index)
                x = F.relu(x)
                x = F.dropout(x, p=0.2, training=self.training)

                x = self.conv2(x, edge_index)
                x = F.relu(x)

                # Global pooling
                x = global_mean_pool(x, batch)

                # Final layer
                x = self.fc(x)
                return x

        # Initialize model
        in_channels = all_data[0].x.size(1)
        model = SimpleGNN(in_channels)

        # Test forward pass with one batch
        print("\nTesting forward pass through GNN...")
        start_time = time.time()

        sample_batch = next(iter(loader))
        output = model(sample_batch.x, sample_batch.edge_index, sample_batch.batch)

        end_time = time.time()
        print(f"Forward pass completed in {end_time - start_time:.3f} seconds")
        print(f"Output shape: {output.shape}")

    except ImportError:
        print("GNN modules not available for testing forward pass")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Load synthetic protein data into PyTorch Geometric")
    parser.add_argument("--data_dir", type=str, default="synthetic_data",
                        help="Directory containing synthetic protein data")
    parser.add_argument("--max_graphs", type=int, default=None,
                        help="Maximum number of graphs to load (default: all)")
    parser.add_argument("--batch_size", type=int, default=32,
                        help="Batch size for testing data loader")
    args = parser.parse_args()

    # Load all graphs
    print("=== Loading Protein Graphs ===")
    _, all_data = load_all_graphs(args.data_dir, args.max_graphs)

    # Test batch loading
    print("\n=== Testing Batch Loading ===")
    test_batch_loading(all_data, args.batch_size)