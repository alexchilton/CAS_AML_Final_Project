import sys
import os
import torch
import numpy as np
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data, Batch
from typing import Dict, List, Any, Tuple, Optional, Union

# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../preprocessing'))

# Import the required classes
from ProteinGraphDataset import ProteinGraphDataset
from GraphBuilder import GraphBuilder
from ProteinFeatureProcessor import ProteinFeatureProcessor


class DatasetHandler:
    """
    Handler for loading and preparing the protein graph dataset.
    Ensures all graph objects have the required attributes and proper batching.
    Supports feature normalization and denormalization at the dataset level.
    """

    def __init__(self, data_root, batch_size=32):
        """
        Initialize the dataset handler.

        Args:
            data_root: Root directory containing protein graph data
            batch_size: Batch size for data loaders
        """
        self.data_root = data_root
        self.batch_size = batch_size
        self.dataset = None
        self.train_loader = None
        self.val_loader = None
        self.test_loader = None

        # Initialize feature processor with default constraints
        self.feature_processor = ProteinFeatureProcessor()

        # Set default feature constraints
        self.set_default_constraints()

    def set_default_constraints(self):
        """
        Set default constraints for the feature processor.
        These constraints ensure that normalized features stay within reasonable bounds.
        """
        # B-factors are typically positive
        self.feature_processor.set_feature_constraints('b_factor', min_val=0.0, max_val=None)

        # Distances between atoms are positive
        self.feature_processor.set_feature_constraints('distance', min_val=0.0, max_val=None)

        # No specific constraints for coordinates by default
        self.feature_processor.set_feature_constraints('coords', min_val=None, max_val=None)

    def set_feature_constraints(self, feature_type: str, min_val: Optional[float] = None,
                                max_val: Optional[float] = None):
        """
        Set constraints for a specific feature type.

        Args:
            feature_type: Type of feature ('coords', 'b_factor', or 'distance')
            min_val: Minimum allowed value (None for no constraint)
            max_val: Maximum allowed value (None for no constraint)
        """
        self.feature_processor.set_feature_constraints(feature_type, min_val, max_val)

    def load_dataset(self, use_global_normalization=True, preprocess_all=False):
        """
        Load the protein graph dataset and set up global normalization if requested.

        Args:
            use_global_normalization: Whether to use dataset-level normalization
            preprocess_all: Whether to preprocess all proteins (needed for global normalization)

        Returns:
            ProteinGraphDataset object
        """
        print("Loading dataset...")

        # Create dataset
        self.dataset = ProteinGraphDataset(
            root=self.data_root,
            #feature_processor=self.feature_processor
        )

        print(f"Dataset loaded with {len(self.dataset)} proteins")

        # If using global normalization, we need to compute statistics first
        if use_global_normalization:
            print("Computing global normalization parameters...")

            # We need to process each graph to collect statistics
            graph_builder = GraphBuilder(self.data_root)
            nanobody_ids = self.dataset.nanobody_ids

            # Process a subset for large datasets
            max_samples = min(len(nanobody_ids), 500)  # Use up to 100 samples
            sample_ids = nanobody_ids[:max_samples]

            print(f"Collecting statistics from {len(sample_ids)} samples...")
            for i, nanobody_id in enumerate(sample_ids):
                if i % 10 == 0:
                    print(f"Processing sample {i+1}/{len(sample_ids)}...")

                # Load graph and collect statistics
                try:
                    graph = graph_builder.load_graph(nanobody_id)
                    self.feature_processor.collect_statistics(graph)
                except Exception as e:
                    print(f"Error processing {nanobody_id}: {e}")

            # Compute global normalization parameters
            global_params = self.feature_processor.compute_global_normalization_params()
            print("Global normalization parameters computed.")

            # Now force preprocessing of all samples with global parameters if requested
            if preprocess_all:
                print("Preprocessing all samples with global normalization...")
                self.dataset.process()

        # Verify and fix dataset
        self._verify_and_fix_dataset()

        # Print dataset statistics
        self._print_statistics()

        return self.dataset

    def set_global_normalization_params(self, params: Dict[str, Any]):
        """
        Set global normalization parameters directly.

        Args:
            params: Dictionary of normalization parameters
        """
        print("Setting global normalization parameters...")
        self.feature_processor.set_global_normalization_params(params)

    def get_global_normalization_params(self) -> Dict[str, Any]:
        """
        Get the current global normalization parameters.

        Returns:
            Dictionary of global normalization parameters
        """
        return self.feature_processor.global_normalization_params

    def _verify_and_fix_dataset(self):
        """
        Verify that all graphs have the required attributes and fix if needed.
        """
        for i in range(len(self.dataset)):
            data = self.dataset[i]

            # Ensure edge_attr exists
            if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                # Create default edge attributes (all ones)
                if hasattr(data, 'edge_index') and data.edge_index.size(1) > 0:
                    data.edge_attr = torch.ones(data.edge_index.size(1), 1, dtype=torch.float)
                else:
                    # If no edges, create empty edge_attr
                    data.edge_attr = torch.empty((0, 1), dtype=torch.float)

            # Ensure node features exist and have proper shape
            if not hasattr(data, 'x') or data.x is None or data.x.nelement() == 0:
                print(f"Warning: Graph {i} has no node features. Creating default features.")
                # Create default node features (all zeros)
                if hasattr(data, 'num_nodes') and data.num_nodes > 0:
                    # Use standard feature dimension (from your ProteinFeatureProcessor)
                    feature_dim = 21 + 7 + 3 + 1 + 7  # amino acids + ss types + coords + b_factor + meiler
                    data.x = torch.zeros(data.num_nodes, feature_dim, dtype=torch.float)
                else:
                    # Empty graph
                    feature_dim = 21 + 7 + 3 + 1 + 7
                    data.x = torch.empty((0, feature_dim), dtype=torch.float)

            # Ensure normalization parameters are present
            if not hasattr(data, 'normalization_params') and self.feature_processor.use_global_normalization:
                data.normalization_params = self.feature_processor.global_normalization_params.copy()

    def _print_statistics(self):
        """
        Print dataset statistics.
        """
        stats = self.dataset.get_statistics()
        print("Dataset statistics:")
        for key, value in stats.items():
            print(f"  {key}: {value}")

    def prepare_dataloaders(self, train_ratio=0.7, val_ratio=0.15):
        """
        Split the dataset and create data loaders.

        Args:
            train_ratio: Ratio of training data
            val_ratio: Ratio of validation data

        Returns:
            Tuple of (train_loader, val_loader, test_loader)
        """
        # Split dataset into train, validation, and test sets
        train_size = int(train_ratio * len(self.dataset))
        val_size = int(val_ratio * len(self.dataset))
        test_size = len(self.dataset) - train_size - val_size

        train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(
            self.dataset, [train_size, val_size, test_size]
        )

        # Create data loaders
        self.train_loader = DataLoader(
            train_dataset,
            batch_size=self.batch_size,
            shuffle=True
        )

        self.val_loader = DataLoader(
            val_dataset,
            batch_size=self.batch_size,
            shuffle=False
        )

        self.test_loader = DataLoader(
            test_dataset,
            batch_size=self.batch_size,
            shuffle=False
        )

        return self.train_loader, self.val_loader, self.test_loader

    def get_input_dim(self):
        """
        Get input dimension from the dataset.

        Returns:
            Input dimension for the model
        """
        sample_data = self.dataset[0]
        return sample_data.x.size(1)

    def denormalize_data(self, data: Union[Data, Batch]) -> Union[Data, Batch]:
        """
        Denormalize the features in a PyTorch Geometric Data or Batch object.

        Args:
            data: PyTorch Geometric Data or Batch object with normalized features

        Returns:
            Data or Batch object with denormalized features
        """
        if not hasattr(data, 'normalization_params'):
            # Try to use global normalization parameters
            if self.feature_processor.use_global_normalization:
                norm_params = self.feature_processor.global_normalization_params
            else:
                raise ValueError("Data object does not have normalization parameters")
        else:
            norm_params = data.normalization_params

        # Convert tensors to numpy for denormalization
        x_np = data.x.cpu().numpy()
        edge_attr_np = data.edge_attr.cpu().numpy() if hasattr(data, 'edge_attr') else None

        # Denormalize node features
        denorm_x_np = self.feature_processor.denormalize_features(x_np, norm_params)

        # Denormalize edge features if present
        if edge_attr_np is not None:
            denorm_edge_attr_np = self.feature_processor.denormalize_edge_features(edge_attr_np, norm_params)
            denorm_edge_attr = torch.tensor(denorm_edge_attr_np, dtype=torch.float)
        else:
            denorm_edge_attr = data.edge_attr

        # Convert back to tensors
        denorm_x = torch.tensor(denorm_x_np, dtype=torch.float)

        # Create new data object with denormalized features
        denorm_data = data.clone()
        denorm_data.x = denorm_x
        if hasattr(data, 'edge_attr'):
            denorm_data.edge_attr = denorm_edge_attr

        return denorm_data

    def save_global_normalization_params(self, path: str):
        """
        Save global normalization parameters to a file.

        Args:
            path: Path to save the parameters
        """
        import json

        params = self.feature_processor.global_normalization_params
        # Convert numpy arrays to lists for JSON serialization
        serializable_params = {}
        for key, value in params.items():
            if isinstance(value, np.ndarray):
                serializable_params[key] = value.tolist()
            else:
                serializable_params[key] = value

        # Create directory if it doesn't exist
        #os.makedirs(os.path.dirname(path), exist_ok=True)

        # Save to file
        with open(path, 'w') as f:
            json.dump(serializable_params, f, indent=4)

        print(f"Global normalization parameters saved to {path}")

    def load_global_normalization_params(self, path: str):
        """
        Load global normalization parameters from a file and apply them.

        Args:
            path: Path to the parameter file
        """
        import json

        with open(path, 'r') as f:
            params = json.load(f)

        # Convert lists back to numpy arrays where needed
        for key in ['coords_mean', 'coords_std']:
            if key in params and isinstance(params[key], list):
                params[key] = np.array(params[key])

        # Set the global normalization parameters
        self.feature_processor.set_global_normalization_params(params)

        print(f"Global normalization parameters loaded from {path}")
        return params