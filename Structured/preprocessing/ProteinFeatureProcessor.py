import networkx as nx
import numpy as np
import torch
from torch_geometric.data import Data
from typing import Dict, List, Tuple, Set, Any, Optional


class ProteinFeatureProcessor:
    """
    A class for processing protein graph features to prepare them for GNN models.
    Works with NetworkX graphs from GraphBuilder to extract and normalize features.
    """

    def __init__(self):
        """Initialize the feature processor with default settings."""
        # Amino acid vocabulary for one-hot encoding
        self.amino_acids = [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL', 'UNK'  # UNK for unknown/non-standard residues
        ]
        self.aa_to_idx = {aa: i for i, aa in enumerate(self.amino_acids)}

        # Secondary structure vocabulary
        self.ss_types = ['H', 'E', 'B', 'G', 'T', 'S', '?']  # From README
        self.ss_to_idx = {ss: i for i, ss in enumerate(self.ss_types)}

    def process_graph_to_pyg(self, graph: nx.Graph) -> Data:
        """
        Process a NetworkX protein graph into a PyTorch Geometric Data object
        with appropriate node and edge features.

        Args:
            graph: NetworkX graph representing a protein structure

        Returns:
            PyTorch Geometric Data object with processed features
        """
        # Extract node features and edge indices
        node_features, node_mapping = self.extract_node_features(graph)
        edge_indices, edge_features = self.extract_edge_features(graph, node_mapping)

        # Convert to PyTorch tensors
        x = torch.tensor(node_features, dtype=torch.float)
        edge_index = torch.tensor(edge_indices, dtype=torch.long)
        edge_attr = torch.tensor(edge_features, dtype=torch.float)

        # Create PyTorch Geometric Data object
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

        return data

    def extract_node_features(self, graph: nx.Graph) -> Tuple[np.ndarray, Dict[str, int]]:
        """
        Extract and normalize node features from a protein graph.

        Args:
            graph: NetworkX graph representing a protein structure

        Returns:
            Tuple of (node_features_array, node_id_to_index_mapping)
        """
        # Create mapping from node IDs to indices
        nodes = list(graph.nodes())
        node_mapping = {node: i for i, node in enumerate(nodes)}

        # Initialize feature array
        num_nodes = len(nodes)
        # Features: amino acid one-hot (21) + secondary structure one-hot (7) +
        #           coordinates (3) + b_factor (1) + meiler features (7)
        feature_dim = len(self.amino_acids) + len(self.ss_types) + 3 + 1 + 7
        node_features = np.zeros((num_nodes, feature_dim))

        # Extract features for each node
        for i, node in enumerate(nodes):
            attrs = graph.nodes[node]

            # Extract residue name and handle non-standard residues
            residue_name = attrs.get('residue_name', 'UNK')
            aa_idx = self.aa_to_idx.get(residue_name, self.aa_to_idx['UNK'])

            # Secondary structure (default to '?' if missing)
            ss = attrs.get('ss', '?')
            ss_idx = self.ss_to_idx.get(ss, self.ss_to_idx['?'])

            # Coordinates (default to origin if missing)
            coords = attrs.get('coords', np.zeros(3))
            if not isinstance(coords, np.ndarray):
                coords = np.array(coords)

            # B-factor (default to 0 if missing)
            b_factor = attrs.get('b_factor', 0.0)

            # Meiler features (default to zeros if missing)
            meiler = attrs.get('meiler', np.zeros(7))
            if not isinstance(meiler, np.ndarray):
                if hasattr(meiler, 'values'):  # Handle pandas Series
                    meiler = meiler.values
                else:
                    meiler = np.array(list(meiler))

            # Create feature vector (one-hot encodings + normalized continuous features)
            feature_idx = 0

            # One-hot encode amino acid
            node_features[i, aa_idx] = 1
            feature_idx += len(self.amino_acids)

            # One-hot encode secondary structure
            node_features[i, feature_idx + ss_idx] = 1
            feature_idx += len(self.ss_types)

            # Add coordinates
            node_features[i, feature_idx:feature_idx+3] = coords[:3]
            feature_idx += 3

            # Add B-factor
            node_features[i, feature_idx] = b_factor
            feature_idx += 1

            # Add Meiler features
            meiler_dim = min(len(meiler), 7)  # Ensure we get at most 7 values
            node_features[i, feature_idx:feature_idx+meiler_dim] = meiler[:meiler_dim]

        # Normalize the continuous features (coordinates, B-factor)
        coord_indices = list(range(len(self.amino_acids) + len(self.ss_types),
                                   len(self.amino_acids) + len(self.ss_types) + 3))
        b_factor_idx = len(self.amino_acids) + len(self.ss_types) + 3
        node_features = self.normalize_features(node_features, coord_indices, b_factor_idx)

        return node_features, node_mapping

    def extract_edge_features(self, graph: nx.Graph, node_mapping: Dict[str, int]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract and normalize edge features from a protein graph.

        Args:
            graph: NetworkX graph representing a protein structure
            node_mapping: Mapping from node IDs to indices

        Returns:
            Tuple of (edge_indices, edge_features)
        """
        # Collect edges and features
        edge_indices = []
        edge_features = []

        for u, v, data in graph.edges(data=True):
            # Add edge indices
            u_idx, v_idx = node_mapping[u], node_mapping[v]
            edge_indices.append([u_idx, v_idx])

            # Get edge features
            features = []

            # Distance feature
            distance = data.get('distance', 0.0)
            features.append(distance)

            # Edge type feature (peptide bond = 1, contact = 0)
            edge_type = 1.0 if 'peptide_bond' in data.get('kind', set()) else 0.0
            features.append(edge_type)

            edge_features.append(features)

        # Create NumPy arrays
        edge_indices = np.array(edge_indices).T  # Transpose to get shape (2, num_edges)
        edge_features = np.array(edge_features)

        # Normalize distance features
        if edge_features.shape[0] > 0:
            distances = edge_features[:, 0]
            if len(distances) > 0:
                mean_dist = np.mean(distances)
                std_dist = np.std(distances)
                if std_dist > 0:
                    edge_features[:, 0] = (distances - mean_dist) / std_dist

        return edge_indices, edge_features

    def normalize_features(self, features: np.ndarray, coord_indices: List[int],
                           b_factor_idx: int = None) -> np.ndarray:
        """
        Normalize continuous features in a feature array.

        Args:
            features: NumPy array of features
            coord_indices: Indices of coordinate features
            b_factor_idx: Index of B-factor feature (if any)

        Returns:
            NumPy array with normalized features
        """
        normalized = features.copy()

        # Normalize 3D coordinates (mean centering)
        if coord_indices:
            coords = normalized[:, coord_indices]
            coords_mean = np.mean(coords, axis=0)
            coords_std = np.std(coords, axis=0)
            coords_std[coords_std == 0] = 1.0  # Avoid division by zero
            normalized[:, coord_indices] = (coords - coords_mean) / coords_std

        # Normalize B-factors if present
        if b_factor_idx is not None:
            b_factors = normalized[:, b_factor_idx]
            b_mean = np.mean(b_factors)
            b_std = np.std(b_factors)
            if b_std > 0:
                normalized[:, b_factor_idx] = (b_factors - b_mean) / b_std

        return normalized