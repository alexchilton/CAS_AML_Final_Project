import networkx as nx
import numpy as np
import torch
from torch_geometric.data import Data
from typing import Dict, List, Tuple, Set, Any, Optional


class ProteinFeatureProcessor:
    """
    A class for processing protein graph features to prepare them for GNN models.
    Works with NetworkX graphs from GraphBuilder to extract and normalize features.
    Includes support for feature constraints and dataset-level normalization.
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

        # Normalization parameters
        self.global_normalization_params = {}
        self.use_global_normalization = False

        # Default feature constraints
        self.feature_constraints = {
            'coords': {'min': None, 'max': None},  # No constraints by default
            'b_factor': {'min': 0.0, 'max': None},  # B-factors typically non-negative
            'distance': {'min': 0.0, 'max': None}   # Distances are non-negative
        }

        # For collecting dataset statistics
        self.dataset_statistics = {
            'coords': {'values': [], 'indices': []},
            'b_factor': {'values': [], 'idx': None},
            'distance': {'values': []}
        }
        self.statistics_collected = False

    def set_feature_constraints(self, feature_type: str, min_val: Optional[float] = None,
                                max_val: Optional[float] = None):
        """
        Set constraints for a feature type.

        Args:
            feature_type: Type of feature ('coords', 'b_factor', or 'distance')
            min_val: Minimum allowed value (None for no constraint)
            max_val: Maximum allowed value (None for no constraint)
        """
        if feature_type in self.feature_constraints:
            self.feature_constraints[feature_type]['min'] = min_val
            self.feature_constraints[feature_type]['max'] = max_val
        else:
            raise ValueError(f"Unknown feature type: {feature_type}")

    def collect_statistics(self, graph: nx.Graph):
        """
        Collect feature statistics from a graph for dataset-level normalization.

        Args:
            graph: NetworkX graph representing a protein structure
        """
        # Extract features without normalization
        nodes = list(graph.nodes())
        num_nodes = len(nodes)
        node_mapping = {node: i for i, node in enumerate(nodes)}

        # Feature dimensions
        feature_dim = len(self.amino_acids) + len(self.ss_types) + 3 + 1 + 7
        node_features = np.zeros((num_nodes, feature_dim))

        # Extract node features
        for i, node in enumerate(nodes):
            attrs = graph.nodes[node]

            # Extract coords
            coords = attrs.get('coords', np.zeros(3))
            if not isinstance(coords, np.ndarray):
                coords = np.array(coords)

            # Extract B-factor
            b_factor = attrs.get('b_factor', 0.0)

            # Calculate feature indices
            coord_start_idx = len(self.amino_acids) + len(self.ss_types)
            b_factor_idx = coord_start_idx + 3

            # Store coordinate indices if not already stored
            if not self.dataset_statistics['coords']['indices']:
                self.dataset_statistics['coords']['indices'] = list(range(coord_start_idx, coord_start_idx + 3))

            # Store B-factor index if not already stored
            if self.dataset_statistics['b_factor']['idx'] is None:
                self.dataset_statistics['b_factor']['idx'] = b_factor_idx

            # Apply constraints
            coords = self._apply_constraint_to_value(coords, 'coords')
            b_factor = self._apply_constraint_to_value(b_factor, 'b_factor')

            # Collect statistics
            self.dataset_statistics['coords']['values'].append(coords)
            self.dataset_statistics['b_factor']['values'].append(b_factor)

        # Extract edge features
        for u, v, data in graph.edges(data=True):
            distance = data.get('distance', 0.0)
            distance = self._apply_constraint_to_value(distance, 'distance')
            self.dataset_statistics['distance']['values'].append(distance)

    def _apply_constraint_to_value(self, value, feature_type):
        """
        Apply constraints to a feature value.

        Args:
            value: Feature value(s) to constrain
            feature_type: Type of feature ('coords', 'b_factor', or 'distance')

        Returns:
            Constrained value(s)
        """
        constraints = self.feature_constraints.get(feature_type, {'min': None, 'max': None})

        if isinstance(value, np.ndarray):
            if constraints['min'] is not None:
                value = np.maximum(value, constraints['min'])
            if constraints['max'] is not None:
                value = np.minimum(value, constraints['max'])
        else:
            if constraints['min'] is not None:
                value = max(value, constraints['min'])
            if constraints['max'] is not None:
                value = min(value, constraints['max'])

        return value

    def compute_global_normalization_params(self):
        """
        Compute global normalization parameters from collected statistics.

        Returns:
            Dictionary of global normalization parameters
        """
        params = {}

        # Process coordinates
        coords_values = np.array(self.dataset_statistics['coords']['values'])
        coords_mean = np.mean(coords_values, axis=0)
        coords_std = np.std(coords_values, axis=0)
        coords_std[coords_std == 0] = 1.0  # Avoid division by zero

        params['coords_mean'] = coords_mean.tolist()
        params['coords_std'] = coords_std.tolist()
        params['coord_indices'] = self.dataset_statistics['coords']['indices']

        # Process B-factors
        b_factor_values = np.array(self.dataset_statistics['b_factor']['values'])
        b_factor_mean = np.mean(b_factor_values)
        b_factor_std = np.std(b_factor_values)
        if b_factor_std == 0:
            b_factor_std = 1.0

        params['b_factor_mean'] = float(b_factor_mean)
        params['b_factor_std'] = float(b_factor_std)
        params['b_factor_idx'] = self.dataset_statistics['b_factor']['idx']

        # Process distances
        distance_values = np.array(self.dataset_statistics['distance']['values'])
        distance_mean = np.mean(distance_values)
        distance_std = np.std(distance_values)
        if distance_std == 0:
            distance_std = 1.0

        params['edge_distance_mean'] = float(distance_mean)
        params['edge_distance_std'] = float(distance_std)

        # Set global parameters
        self.global_normalization_params = params
        self.statistics_collected = True
        self.use_global_normalization = True

        # Reset statistics storage to free memory
        self.dataset_statistics = {
            'coords': {'values': [], 'indices': self.dataset_statistics['coords']['indices']},
            'b_factor': {'values': [], 'idx': self.dataset_statistics['b_factor']['idx']},
            'distance': {'values': []}
        }

        return params

    def set_global_normalization_params(self, params: Dict[str, Any]):
        """
        Set global normalization parameters directly.

        Args:
            params: Dictionary of normalization parameters
        """
        self.global_normalization_params = params
        self.use_global_normalization = True
        self.statistics_collected = True

    def process_graph_to_pyg(self, graph: nx.Graph) -> Data:
        """
        Process a NetworkX protein graph into a PyTorch Geometric Data object
        with appropriate node and edge features.

        Args:
            graph: NetworkX graph representing a protein structure

        Returns:
            PyTorch Geometric Data object with processed features
        """
        # If in statistics collection mode and global params not computed yet
        if not self.statistics_collected and not self.use_global_normalization:
            self.collect_statistics(graph)

        # Extract node features and edge indices
        node_features, node_mapping = self.extract_node_features(graph)
        edge_indices, edge_features = self.extract_edge_features(graph, node_mapping)

        # Convert to PyTorch tensors
        x = torch.tensor(node_features, dtype=torch.float)
        edge_index = torch.tensor(edge_indices, dtype=torch.long)
        edge_attr = torch.tensor(edge_features, dtype=torch.float)

        # Create PyTorch Geometric Data object
        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)

        # Store normalization parameters in the data object
        if self.use_global_normalization:
            data.normalization_params = self.global_normalization_params

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

        # Get coordinate indices and B-factor index
        coord_indices = list(range(len(self.amino_acids) + len(self.ss_types),
                                   len(self.amino_acids) + len(self.ss_types) + 3))
        b_factor_idx = len(self.amino_acids) + len(self.ss_types) + 3

        # Apply constraints before normalization
        node_features = self._apply_constraints(node_features,
                                                coords_indices=coord_indices,
                                                b_factor_idx=b_factor_idx)

        # Normalize with parameter storage
        if self.use_global_normalization:
            node_features = self.apply_global_normalization(node_features, "node")
        else:
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

        # Apply constraints to distance features
        if edge_features.shape[0] > 0:
            distances = edge_features[:, 0]
            constraints = self.feature_constraints['distance']
            if constraints['min'] is not None:
                distances = np.maximum(distances, constraints['min'])
            if constraints['max'] is not None:
                distances = np.minimum(distances, constraints['max'])
            edge_features[:, 0] = distances

        # Normalize distance features
        if edge_features.shape[0] > 0:
            if self.use_global_normalization:
                edge_features = self.apply_global_normalization(edge_features, "edge")
            else:
                distances = edge_features[:, 0]
                if len(distances) > 0:
                    mean_dist = np.mean(distances)
                    std_dist = np.std(distances)
                    if std_dist > 0:
                        edge_features[:, 0] = (distances - mean_dist) / std_dist

                    # If not using global normalization, store per-graph parameters
                    if not self.use_global_normalization:
                        self.global_normalization_params['edge_distance_mean'] = float(mean_dist)
                        self.global_normalization_params['edge_distance_std'] = float(std_dist)

        return edge_indices, edge_features

    def apply_global_normalization(self, features: np.ndarray, feature_type: str) -> np.ndarray:
        """
        Apply global normalization parameters to features.

        Args:
            features: Features to normalize
            feature_type: Type of features ('node' or 'edge')

        Returns:
            Normalized features
        """
        normalized = features.copy()

        if feature_type == "node":
            # Normalize coordinates
            if 'coord_indices' in self.global_normalization_params:
                indices = self.global_normalization_params['coord_indices']
                mean = np.array(self.global_normalization_params['coords_mean'])
                std = np.array(self.global_normalization_params['coords_std'])
                normalized[:, indices] = (features[:, indices] - mean) / std

            # Normalize B-factor
            if 'b_factor_idx' in self.global_normalization_params:
                idx = self.global_normalization_params['b_factor_idx']
                mean = self.global_normalization_params['b_factor_mean']
                std = self.global_normalization_params['b_factor_std']
                normalized[:, idx] = (features[:, idx] - mean) / std

        elif feature_type == "edge" and features.shape[0] > 0:
            # Normalize distance
            if 'edge_distance_mean' in self.global_normalization_params:
                mean = self.global_normalization_params['edge_distance_mean']
                std = self.global_normalization_params['edge_distance_std']
                normalized[:, 0] = (features[:, 0] - mean) / std

        return normalized

    def normalize_features(self, features: np.ndarray, coord_indices: List[int],
                           b_factor_idx: int = None) -> np.ndarray:
        """
        Normalize continuous features in a feature array and store normalization parameters.
        Used when global normalization is not enabled.

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

            # Store normalization parameters
            self.global_normalization_params['coords_mean'] = coords_mean.tolist()
            self.global_normalization_params['coords_std'] = coords_std.tolist()
            self.global_normalization_params['coord_indices'] = coord_indices

        # Normalize B-factors if present
        if b_factor_idx is not None:
            b_factors = normalized[:, b_factor_idx]
            b_mean = np.mean(b_factors)
            b_std = np.std(b_factors)
            if b_std > 0:
                normalized[:, b_factor_idx] = (b_factors - b_mean) / b_std

            # Store normalization parameters
            self.global_normalization_params['b_factor_mean'] = float(b_mean)
            self.global_normalization_params['b_factor_std'] = float(b_std)
            self.global_normalization_params['b_factor_idx'] = b_factor_idx

        return normalized

    def denormalize_features(self, normalized_features: np.ndarray,
                             normalization_params: Dict[str, Any]) -> np.ndarray:
        """
        Denormalize features using stored normalization parameters.

        Args:
            normalized_features: NumPy array of normalized features
            normalization_params: Dictionary of normalization parameters

        Returns:
            NumPy array with denormalized features
        """
        denormalized = normalized_features.copy()

        # Denormalize coordinates
        if 'coord_indices' in normalization_params:
            indices = normalization_params['coord_indices']
            mean = np.array(normalization_params['coords_mean'])
            std = np.array(normalization_params['coords_std'])
            denormalized[:, indices] = normalized_features[:, indices] * std + mean

        # Denormalize B-factor
        if 'b_factor_idx' in normalization_params:
            idx = normalization_params['b_factor_idx']
            mean = normalization_params['b_factor_mean']
            std = normalization_params['b_factor_std']
            denormalized[:, idx] = normalized_features[:, idx] * std + mean

        return denormalized

    def denormalize_edge_features(self, normalized_edge_features: np.ndarray,
                                  normalization_params: Dict[str, Any]) -> np.ndarray:
        """
        Denormalize edge features using stored normalization parameters.

        Args:
            normalized_edge_features: NumPy array of normalized edge features
            normalization_params: Dictionary of normalization parameters

        Returns:
            NumPy array with denormalized edge features
        """
        denormalized = normalized_edge_features.copy()

        # Denormalize distance features
        if 'edge_distance_mean' in normalization_params and 'edge_distance_std' in normalization_params:
            mean = normalization_params['edge_distance_mean']
            std = normalization_params['edge_distance_std']
            denormalized[:, 0] = normalized_edge_features[:, 0] * std + mean

        return denormalized

    def _apply_constraints(self, features: np.ndarray, coords_indices: List[int],
                           b_factor_idx: int = None) -> np.ndarray:
        """
        Apply constraints to features before normalization.

        Args:
            features: NumPy array of features
            coords_indices: Indices of coordinate features
            b_factor_idx: Index of B-factor feature (if any)

        Returns:
            NumPy array with constrained features
        """
        constrained = features.copy()

        # Apply constraints to coordinates
        if coords_indices:
            constraints = self.feature_constraints['coords']
            coords = constrained[:, coords_indices]

            if constraints['min'] is not None:
                coords = np.maximum(coords, constraints['min'])
            if constraints['max'] is not None:
                coords = np.minimum(coords, constraints['max'])

            constrained[:, coords_indices] = coords

        # Apply constraints to B-factor
        if b_factor_idx is not None:
            constraints = self.feature_constraints['b_factor']
            b_factors = constrained[:, b_factor_idx]

            if constraints['min'] is not None:
                b_factors = np.maximum(b_factors, constraints['min'])
            if constraints['max'] is not None:
                b_factors = np.minimum(b_factors, constraints['max'])

            constrained[:, b_factor_idx] = b_factors

        return constrained