import os
import torch
import numpy as np
from torch_geometric.data import Dataset, Data
from typing import List, Dict, Any, Optional, Tuple, Union, Callable
import networkx as nx
import pickle
import json

from GraphBuilder import GraphBuilder
from ProteinFeatureProcessor import ProteinFeatureProcessor


class ProteinGraphDataset(Dataset):
    """
    PyTorch Geometric Dataset for protein graph structures.
    Loads and processes protein graphs from the specified directory.
    """

    def __init__(self, root: str, transform: Optional[Callable] = None,
                 pre_transform: Optional[Callable] = None,
                 use_cache: bool = True,
                 selected_nanobodies: Optional[List[str]] = None):
        """
        Initialize the protein graph dataset.

        Args:
            root: Root directory containing protein graph data
            transform: Transform to apply to each data object
            pre_transform: Transform to apply to all data before saving
            use_cache: Whether to cache processed data to disk
            selected_nanobodies: Optional list of specific nanobody IDs to include
        """
        self.root_dir = root
        self.use_cache = use_cache
        self.selected_nanobodies = selected_nanobodies

        # Initialize GraphBuilder and FeatureProcessor
        self.graph_builder = GraphBuilder(root)
        self.feature_processor = ProteinFeatureProcessor()

        # Get list of available nanobody IDs
        self.nanobody_ids = self._filter_nanobodies(self.graph_builder.get_nanobody_ids())

        # Create processed directory if it doesn't exist and caching is enabled
        if self.use_cache and not os.path.exists(self.processed_dir):
            os.makedirs(self.processed_dir)

        # Call parent class initialization at the end
        # This lets PyTorch Geometric set up directory structure
        super().__init__(root, transform, pre_transform)

    @property
    def processed_dir(self) -> str:
        """
        Override the processed_dir property to use a directory in the root.

        Returns:
            Path to the processed directory
        """
        return os.path.join(self.root_dir, 'processed')

    def _filter_nanobodies(self, nanobody_ids: List[str]) -> List[str]:
        """
        Filter nanobody IDs based on selection criteria.

        Args:
            nanobody_ids: List of all available nanobody IDs

        Returns:
            Filtered list of nanobody IDs
        """
        if self.selected_nanobodies is not None:
            # Filter to include only selected nanobodies
            return [n_id for n_id in nanobody_ids if n_id in self.selected_nanobodies]
        return nanobody_ids

    @property
    def raw_file_names(self) -> List[str]:
        """
        Get list of raw file names.

        Returns:
            List of raw file paths (empty list as we handle raw files ourselves)
        """
        return []

    @property
    def processed_file_names(self) -> List[str]:
        """
        Get list of processed file names.

        Returns:
            List of processed file paths
        """
        return [f"{nanobody_id}.pt" for nanobody_id in self.nanobody_ids]

    def len(self) -> int:
        """
        Get number of graphs in the dataset.

        Returns:
            Number of graphs
        """
        return len(self.nanobody_ids)

    def get(self, idx: int) -> Data:
        """
        Get a specific graph by index.

        Args:
            idx: Index of the graph to retrieve

        Returns:
            PyTorch Geometric Data object
        """
        nanobody_id = self.nanobody_ids[idx]

        # Check if processed file exists
        processed_file = os.path.join(self.processed_dir, f"{nanobody_id}.pt")

        if self.use_cache and os.path.exists(processed_file):
            # Load processed data from cache
            data = torch.load(processed_file)
        else:
            # Process data from scratch
            data = self._process_nanobody(nanobody_id)

            # Save processed data to cache
            if self.use_cache:
                torch.save(data, processed_file)

        return data

    def _process_nanobody(self, nanobody_id: str) -> Data:
        """
        Process a single nanobody into a PyTorch Geometric Data object.

        Args:
            nanobody_id: ID of the nanobody to process

        Returns:
            PyTorch Geometric Data object
        """
        try:
            # Load graph using GraphBuilder
            graph = self.graph_builder.load_graph(nanobody_id)

            # Process graph using FeatureProcessor
            data = self.feature_processor.process_graph_to_pyg(graph)

            # Add nanobody ID as attribute
            data.nanobody_id = nanobody_id

            return data

        except Exception as e:
            print(f"Error processing nanobody {nanobody_id}: {str(e)}")
            # Return empty Data object in case of error
            return Data(x=torch.tensor([], dtype=torch.float),
                        edge_index=torch.tensor([[],[]], dtype=torch.long),
                        nanobody_id=nanobody_id)

    def process(self):
        """
        Process all nanobodies and save to disk for future use.
        """
        if not self.use_cache:
            return

        # Process all nanobodies and save to disk
        for i, nanobody_id in enumerate(self.nanobody_ids):
            print(f"Processing {i+1}/{len(self.nanobody_ids)}: {nanobody_id}")

            processed_file = os.path.join(self.processed_dir, f"{nanobody_id}.pt")

            # Skip if already processed
            if os.path.exists(processed_file):
                continue

            # Process nanobody
            data = self._process_nanobody(nanobody_id)

            # Save processed data
            torch.save(data, processed_file)

    def get_statistics(self) -> Dict[str, Any]:
        """
        Compute statistics about the dataset.

        Returns:
            Dictionary of dataset statistics
        """
        num_nodes = []
        num_edges = []

        # Sample up to 100 graphs for statistics
        sample_size = min(100, len(self))
        indices = np.random.choice(len(self), sample_size, replace=False)

        for idx in indices:
            data = self.get(idx)
            num_nodes.append(data.x.shape[0])
            num_edges.append(data.edge_index.shape[1])

        stats = {
            "total_graphs": len(self),
            "avg_nodes": np.mean(num_nodes),
            "std_nodes": np.std(num_nodes),
            "min_nodes": np.min(num_nodes),
            "max_nodes": np.max(num_nodes),
            "avg_edges": np.mean(num_edges),
            "std_edges": np.std(num_edges),
            "min_edges": np.min(num_edges),
            "max_edges": np.max(num_edges),
        }

        return stats