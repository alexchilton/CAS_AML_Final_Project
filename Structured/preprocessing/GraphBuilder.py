import os
import pickle
import networkx as nx
import json
from typing import List, Dict, Any, Optional, Set, Tuple
import torch
from torch_geometric.data import Data


class GraphBuilder:
    """
    A class for loading and processing protein structure data from NetworkX graphs
    and associated pickle files into PyTorch Geometric format.
    """

    def __init__(self, base_dir: str):
        """
        Initialize the GraphBuilder with the base directory containing protein data.

        Args:
            base_dir: Path to the directory containing protein data subdirectories
        """
        self.base_dir = base_dir
        self.nanobody_dirs = self._find_nanobody_dirs()

    def _find_nanobody_dirs(self) -> List[str]:
        """
        Find all nanobody directories in the base directory.

        Returns:
            List of directory names for nanobodies
        """
        if not os.path.exists(self.base_dir):
            raise FileNotFoundError(f"Base directory {self.base_dir} not found")

        return [d for d in os.listdir(self.base_dir)
                if os.path.isdir(os.path.join(self.base_dir, d))]

    def get_nanobody_ids(self) -> List[str]:
        """
        Get the list of available nanobody IDs.

        Returns:
            List of nanobody IDs
        """
        return self.nanobody_dirs

    def load_graph(self, nanobody_id: str) -> nx.Graph:
        """
        Load the NetworkX graph for a specific nanobody.

        Args:
            nanobody_id: ID of the nanobody to load

        Returns:
            NetworkX graph object
        """
        if nanobody_id not in self.nanobody_dirs:
            raise ValueError(f"Nanobody ID {nanobody_id} not found")

        graph_path = os.path.join(self.base_dir, nanobody_id, f"{nanobody_id}_graph.pkl")

        if not os.path.exists(graph_path):
            raise FileNotFoundError(f"Graph file not found at {graph_path}")

        with open(graph_path, 'rb') as f:
            graph = pickle.load(f)

        return graph

    def load_data_dict(self, nanobody_id: str) -> Dict[str, Any]:
        """
        Load the main data dictionary for a specific nanobody.

        Args:
            nanobody_id: ID of the nanobody to load

        Returns:
            Dictionary containing protein data
        """
        if nanobody_id not in self.nanobody_dirs:
            raise ValueError(f"Nanobody ID {nanobody_id} not found")

        data_path = os.path.join(self.base_dir, nanobody_id, f"{nanobody_id}_data.pkl")

        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data file not found at {data_path}")

        with open(data_path, 'rb') as f:
            data_dict = pickle.load(f)

        return data_dict

    def list_chain_files(self, nanobody_id: str) -> Dict[str, List[str]]:
        """
        Get a list of available chain-specific files for a nanobody.

        Args:
            nanobody_id: ID of the nanobody to check

        Returns:
            Dictionary mapping chain IDs to lists of available files
        """
        if nanobody_id not in self.nanobody_dirs:
            raise ValueError(f"Nanobody ID {nanobody_id} not found")

        nanobody_dir = os.path.join(self.base_dir, nanobody_id)
        files = os.listdir(nanobody_dir)

        # Group files by chain
        chain_files = {}
        for filename in files:
            # Check if the filename starts with nanobody_id and contains additional parts
            if filename.startswith(f"{nanobody_id}_") and '_' in filename[len(nanobody_id)+1:]:
                # Extract chain_id and file_type
                remaining = filename[len(nanobody_id)+1:]
                parts = remaining.split('_')
                if len(parts) >= 1 and parts[0] not in ['graph', 'data']:
                    chain_id = parts[0]
                    if chain_id not in chain_files:
                        chain_files[chain_id] = []
                    chain_files[chain_id].append(filename)

        return chain_files

    def load_chain_file(self, nanobody_id: str, chain_id: str, file_type: str) -> Any:
        """
        Load a specific chain file for a nanobody.

        Args:
            nanobody_id: ID of the nanobody
            chain_id: Chain identifier
            file_type: Type of file to load (e.g., 'atoms', 'backbone', 'ss')

        Returns:
            Content of the requested file
        """
        if nanobody_id not in self.nanobody_dirs:
            raise ValueError(f"Nanobody ID {nanobody_id} not found")

        file_path = os.path.join(self.base_dir, nanobody_id, f"{nanobody_id}_{chain_id}_{file_type}.pkl")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Chain file not found at {file_path}")

        with open(file_path, 'rb') as f:
            data = pickle.load(f)

        return data

    def get_available_file_types(self, nanobody_id: str, chain_id: str) -> List[str]:
        """
        Get a list of available file types for a specific chain.

        Args:
            nanobody_id: ID of the nanobody
            chain_id: Chain identifier

        Returns:
            List of available file types
        """
        chain_files = self.list_chain_files(nanobody_id)

        if chain_id not in chain_files:
            raise ValueError(f"Chain ID {chain_id} not found for nanobody {nanobody_id}")

        pdb_id = nanobody_id.split('_')[0]
        file_types = []

        for filename in chain_files[chain_id]:
            # Extract file type from filename
            file_type = filename.replace(f"{pdb_id}_{chain_id}_", "").replace(".pkl", "")
            if file_type:
                file_types.append(file_type)

        return file_types

    def print_graph_summary(self, graph: nx.Graph) -> None:
        """
        Print a summary of the graph properties.

        Args:
            graph: NetworkX graph object
        """
        print(f"Graph Summary:")
        print(f"  Number of nodes: {graph.number_of_nodes()}")
        print(f"  Number of edges: {graph.number_of_edges()}")

        if graph.nodes:
            # Show sample node attributes
            sample_node = list(graph.nodes)[0]
            print(f"\nSample Node: {sample_node}")
            print("  Attributes:")
            for key, value in graph.nodes[sample_node].items():
                print(f"    {key}: {type(value).__name__}")

        if graph.edges:
            # Show sample edge attributes
            sample_edge = list(graph.edges)[0]
            print(f"\nSample Edge: {sample_edge}")
            print("  Attributes:")
            for key, value in graph.edges[sample_edge].items():
                print(f"    {key}: {type(value).__name__}")

    def convert_to_pytorch_geometric(self, graph: nx.Graph) -> Data:
        """
        Convert a NetworkX graph to PyTorch Geometric format.

        Args:
            graph: NetworkX graph object

        Returns:
            PyTorch Geometric Data object
        """
        # Node features: use coordinates as default features
        nodes = list(graph.nodes())
        node_mapping = {node: i for i, node in enumerate(nodes)}

        # Use CA coordinates as node features if available
        node_features = []
        for node in nodes:
            if 'coords' in graph.nodes[node]:
                node_features.append(graph.nodes[node]['coords'])
            elif 'CA' in graph.nodes[node]:
                node_features.append(graph.nodes[node]['CA'])
            else:
                # Default to zeros if no coordinates found
                node_features.append([0.0, 0.0, 0.0])

        # Edge indices
        edge_indices = []
        edge_attrs = []

        for u, v, data in graph.edges(data=True):
            edge_indices.append([node_mapping[u], node_mapping[v]])
            # Use distance as edge attribute if available
            if 'distance' in data:
                edge_attrs.append([data['distance']])
            else:
                edge_attrs.append([1.0])  # Default edge attribute

        # Convert to PyTorch tensors
        import numpy as np
        x = torch.tensor(np.array(node_features), dtype=torch.float)
        edge_index = torch.tensor(np.array(edge_indices), dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(np.array(edge_attrs), dtype=torch.float)

        # Create Data object
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

        return data