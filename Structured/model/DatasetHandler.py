import sys
import os
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data, Batch

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

    def load_dataset(self):
        """
        Load the protein graph dataset and verify its structure.
        """
        print("Loading dataset...")
        self.dataset = ProteinGraphDataset(root=self.data_root)
        print(f"Dataset loaded with {len(self.dataset)} proteins")

        # Verify and fix dataset
        self._verify_and_fix_dataset()

        # Print dataset statistics
        self._print_statistics()

        return self.dataset

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

            # Update the dataset
            #self.dataset[i] = data

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
        """
        # Split dataset into train, validation, and test sets
        train_size = int(train_ratio * len(self.dataset))
        val_size = int(val_ratio * len(self.dataset))
        test_size = len(self.dataset) - train_size - val_size

        train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(
            self.dataset, [train_size, val_size, test_size]
        )

        # Custom collate function to handle missing attributes
        def collate_fn(batch):
            # Check if all items have the required attributes
            for item in batch:
                if not hasattr(item, 'edge_attr') or item.edge_attr is None:
                    item.edge_attr = torch.ones(item.edge_index.size(1), 1, dtype=torch.float)

            # Use Batch.from_data_list for batching
            return Batch.from_data_list(batch)

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