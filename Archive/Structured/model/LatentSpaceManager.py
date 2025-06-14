import torch
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from typing import List, Dict, Any, Optional
import os
from torch_geometric.loader import DataLoader as PyGDataLoader

from Structured.model  import GATVAE


class LatentSpaceManager:
    """
    Class for analyzing and visualizing the latent space of the GAT-VAE model.
    """

    def __init__(self, model: GATVAE, device='cpu', save_dir='./latent_plots'):
        """
        Initialize the latent space manager.

        Args:
            model: Trained GATVAE model
            device: Device to run computations on
            save_dir: Directory to save visualizations
        """
        self.model = model
        self.device = device
        self.save_dir = save_dir

        # Create directory if it doesn't exist
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Move model to device and set to evaluation mode
        self.model.to(device)
        self.model.eval()

    def encode_dataset(self, dataset_or_loader):
        """
        Encode an entire dataset into the latent space.

        Args:
            dataset_or_loader: Dataset or DataLoader containing the dataset

        Returns:
            latent_vectors: List of latent vectors
            nanobody_ids: List of corresponding nanobody IDs
        """
        latent_vectors = []
        nanobody_ids = []

        # Check if we got a dataset or a dataloader
        if not hasattr(dataset_or_loader, '__iter__') or isinstance(dataset_or_loader, list):
            # This is likely a dataset, create a proper PyG DataLoader
            dataloader = PyGDataLoader(dataset_or_loader, batch_size=8, shuffle=False)
        else:
            # This is already a dataloader
            dataloader = dataset_or_loader

        with torch.no_grad():
            for batch in dataloader:
                # Move data to device
                batch = batch.to(self.device)

                # Encode data
                z = self.model.encode(batch)

                # Get nanobody IDs (handle batched data)
                if hasattr(batch, 'nanobody_id'):
                    if isinstance(batch.nanobody_id, list):
                        # List of IDs (typical for batched data)
                        batch_ids = batch.nanobody_id
                    else:
                        # Single ID
                        batch_ids = [batch.nanobody_id]
                else:
                    # If no nanobody_id attribute, use indices
                    batch_size = batch.num_graphs
                    batch_ids = [f"unknown_{i}" for i in range(batch_size)]

                # For each graph in the batch
                if batch.batch is not None:
                    # Get number of graphs in batch
                    num_graphs = batch.batch.max().item() + 1

                    # Split the latent vectors by graph
                    for i in range(num_graphs):
                        latent_vectors.append(z[i].cpu().detach().numpy())
                else:
                    # Single graph
                    latent_vectors.append(z.cpu().detach().numpy())

                # Add IDs
                nanobody_ids.extend(batch_ids)

        return latent_vectors, nanobody_ids

    def visualize_2d(self, latent_vectors, labels=None, method='pca', show=True, save=True):
        """
        Visualize latent space in 2D using PCA or t-SNE.

        Args:
            latent_vectors: List of latent vectors
            labels: Optional list of labels for coloring points
            method: Dimensionality reduction method ('pca' or 'tsne')
            show: Whether to display the plot
            save: Whether to save the plot to disk
        """
        # Stack latent vectors and convert to numpy
        if isinstance(latent_vectors, list):
            latent_vectors = np.vstack(latent_vectors)

        # Reduce dimensionality
        if method.lower() == 'pca':
            reducer = PCA(n_components=2)
            reduced_vecs = reducer.fit_transform(latent_vectors)
            title = 'PCA Visualization of Latent Space'
        elif method.lower() == 'tsne':
            reducer = TSNE(n_components=2, random_state=42)
            reduced_vecs = reducer.fit_transform(latent_vectors)
            title = 't-SNE Visualization of Latent Space'
        else:
            raise ValueError(f"Unknown method: {method}. Choose 'pca' or 'tsne'.")

        # Create plot
        plt.figure(figsize=(10, 8))

        if labels is not None:
            # Plot with colors based on labels
            scatter = plt.scatter(reduced_vecs[:, 0], reduced_vecs[:, 1],
                                  c=labels, cmap='viridis', alpha=0.7)
            plt.colorbar(scatter, label='Label')
        else:
            # Plot without labels
            plt.scatter(reduced_vecs[:, 0], reduced_vecs[:, 1], alpha=0.7)

        plt.title(title)
        plt.xlabel('Component 1')
        plt.ylabel('Component 2')
        plt.grid(True, alpha=0.3)

        if save:
            plt.savefig(f"{self.save_dir}/latent_space_{method}.png", dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

        return reduced_vecs

    def interpolate(self, data1, data2, steps=10):
        """
        Interpolate between two protein structures in latent space.

        Args:
            data1: First protein graph
            data2: Second protein graph
            steps: Number of interpolation steps

        Returns:
            interpolated_graphs: List of reconstructed graphs
        """
        self.model.eval()

        with torch.no_grad():
            # Move data to device
            data1 = data1.to(self.device)
            data2 = data2.to(self.device)

            # Encode proteins
            z1 = self.model.encode(data1)
            z2 = self.model.encode(data2)

            # Interpolate in latent space
            interpolated_zs = []
            for alpha in np.linspace(0, 1, steps):
                z_interp = (1 - alpha) * z1 + alpha * z2
                interpolated_zs.append(z_interp)

            # Decode interpolated vectors
            interpolated_graphs = []
            for z in interpolated_zs:
                # Use the graph structure from data1
                x_reconstructed = self.model.decode(z, data1)
                # Create a new Data object with reconstructed features
                reconstructed_graph = type(data1)(
                    x=x_reconstructed,
                    edge_index=data1.edge_index,
                    edge_attr=data1.edge_attr
                )
                interpolated_graphs.append(reconstructed_graph)

        return interpolated_graphs

    def sample_from_latent(self, num_samples=5, template_data=None):
        """
        Generate new protein structures by sampling from the latent space.

        Args:
            num_samples: Number of samples to generate
            template_data: Template graph structure to use for reconstruction

        Returns:
            sampled_graphs: List of sampled graphs
        """
        if template_data is None:
            raise ValueError("Template data is required for reconstruction")

        self.model.eval()
        template_data = template_data.to(self.device)

        with torch.no_grad():
            # Sample from normal distribution
            z_samples = self.model.sample(num_samples, device=self.device)

            # Decode samples
            sampled_graphs = []
            for z in z_samples:
                x_reconstructed = self.model.decode(z.unsqueeze(0), template_data)
                # Create a new Data object with reconstructed features
                sampled_graph = type(template_data)(
                    x=x_reconstructed,
                    edge_index=template_data.edge_index,
                    edge_attr=template_data.edge_attr
                )
                sampled_graphs.append(sampled_graph)

        return sampled_graphs