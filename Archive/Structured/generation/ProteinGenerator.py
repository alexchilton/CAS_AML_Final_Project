import sys

import torch
import torch.nn.functional as F
import numpy as np
import os
from typing import List, Dict, Any, Optional, Tuple
from torch_geometric.data import Data
# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../model'))
from GATVAE import GATVAE


class ProteinGenerator:
    """
    Class for generating protein structures from the latent space.
    Handles sampling, interpolation, and conditional generation.
    """

    def __init__(self, model: GATVAE, device='cpu'):
        """
        Initialize the protein generator.

        Args:
            model: Trained GATVAE model
            device: Device to run computations on
        """
        self.model = model
        self.device = device

        # Move model to device and set to evaluation mode
        self.model.to(device)
        self.model.eval()

    def sample_random(self, num_samples: int, template_data: Data) -> List[Data]:
        """
        Generate random protein structures by sampling from the latent space.

        Args:
            num_samples: Number of samples to generate
            template_data: Template graph structure to use for reconstruction

        Returns:
            List of sampled protein graphs
        """
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
                sampled_graph = Data(
                    x=x_reconstructed,
                    edge_index=template_data.edge_index,
                    edge_attr=template_data.edge_attr
                )
                sampled_graphs.append(sampled_graph)

        return sampled_graphs

    def interpolate(self, data1: Data, data2: Data, steps: int = 10) -> List[Data]:
        """
        Interpolate between two protein structures in latent space.

        Args:
            data1: First protein graph
            data2: Second protein graph
            steps: Number of interpolation steps

        Returns:
            List of interpolated protein graphs
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
                reconstructed_graph = Data(
                    x=x_reconstructed,
                    edge_index=data1.edge_index,
                    edge_attr=data1.edge_attr
                )
                interpolated_graphs.append(reconstructed_graph)

        return interpolated_graphs

    def reconstruct(self, data: Data) -> Data:
        """
        Reconstruct a protein graph through the autoencoder.

        Args:
            data: Input protein graph

        Returns:
            Reconstructed protein graph
        """
        self.model.eval()
        data = data.to(self.device)

        with torch.no_grad():
            # Encode
            z = self.model.encode(data)

            # Decode
            x_reconstructed = self.model.decode(z, data)

            # Create new Data object with reconstructed features
            reconstructed = Data(
                x=x_reconstructed,
                edge_index=data.edge_index,
                edge_attr=data.edge_attr
            )

        return reconstructed

    def conditional_generation(self, condition_vector: torch.Tensor, template_data: Data,
                               num_samples: int = 5, std_scale: float = 0.5) -> List[Data]:
        """
        Generate proteins conditioned on a specific latent vector with variations.

        Args:
            condition_vector: Base latent vector to condition on
            template_data: Template graph structure to use for reconstruction
            num_samples: Number of variations to generate
            std_scale: Scale of the random perturbation (0 = no randomness, 1 = full random)

        Returns:
            List of conditionally generated protein graphs
        """
        self.model.eval()
        template_data = template_data.to(self.device)
        condition_vector = condition_vector.to(self.device)

        with torch.no_grad():
            # Generate variations around the condition vector
            sampled_graphs = []

            for _ in range(num_samples):
                # Add controlled random noise to the condition vector
                noise = torch.randn_like(condition_vector) * std_scale
                z = condition_vector + noise

                # Decode
                x_reconstructed = self.model.decode(z.unsqueeze(0), template_data)

                # Create a new Data object with reconstructed features
                sampled_graph = Data(
                    x=x_reconstructed,
                    edge_index=template_data.edge_index,
                    edge_attr=template_data.edge_attr
                )
                sampled_graphs.append(sampled_graph)

        return sampled_graphs

    def average_proteins(self, protein_data_list: List[Data]) -> Data:
        """
        Create a new protein by averaging multiple proteins in latent space.

        Args:
            protein_data_list: List of protein graphs to average

        Returns:
            Averaged protein graph
        """
        self.model.eval()

        if not protein_data_list:
            raise ValueError("Empty protein list provided")

        with torch.no_grad():
            # Use the first protein as template for structure
            template_data = protein_data_list[0].to(self.device)

            # Encode all proteins to get latent vectors
            latent_vectors = []
            for protein in protein_data_list:
                protein = protein.to(self.device)
                z = self.model.encode(protein)
                latent_vectors.append(z)

            # Average the latent vectors
            avg_latent = torch.mean(torch.stack(latent_vectors), dim=0)

            # Decode the average latent vector
            x_reconstructed = self.model.decode(avg_latent.unsqueeze(0), template_data)

            # Create a new Data object with reconstructed features
            averaged_protein = Data(
                x=x_reconstructed,
                edge_index=template_data.edge_index,
                edge_attr=template_data.edge_attr
            )

        return averaged_protein