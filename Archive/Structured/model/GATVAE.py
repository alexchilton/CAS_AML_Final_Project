import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data

from Structured.model import GATEncoder
from Structured.model import GATDecoder


class GATVAE(nn.Module):
    """
    Improved Graph Attention Variational Autoencoder for protein graphs.
    Combines the GATEncoder and ImprovedGATDecoder into a single model.
    """

    def __init__(self, input_dim, hidden_dim=64, latent_dim=32, num_heads=4, dropout=0.1):
        """
        Initialize the GAT-based VAE.

        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            latent_dim: Dimension of the latent space
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
        """
        super(GATVAE, self).__init__()

        # Initialize encoder
        self.encoder = GATEncoder(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            latent_dim=latent_dim,
            num_heads=num_heads,
            dropout=dropout
        )

        # Initialize improved decoder
        self.decoder = GATDecoder(
            latent_dim=latent_dim,
            hidden_dim=hidden_dim,
            output_dim=input_dim,
            num_heads=num_heads,
            dropout=dropout
        )

        self.latent_dim = latent_dim

    def forward(self, data):
        """
        Forward pass through the VAE.

        Args:
            data: PyTorch Geometric Data object containing the graph

        Returns:
            x_reconstructed: Reconstructed node features
            mu: Mean of the latent distribution
            logvar: Log variance of the latent distribution
        """
        # Encode
        mu, logvar = self.encoder(data)

        # Sample from latent distribution
        z = self.encoder.reparameterize(mu, logvar)

        # Decode
        x_reconstructed = self.decoder(z, data)

        return x_reconstructed, mu, logvar

    def encode(self, data):
        """
        Encode data to latent representation.

        Args:
            data: PyTorch Geometric Data object

        Returns:
            z: Latent vector
        """
        return self.encoder.encode(data)

    def decode(self, z, data):
        """
        Decode latent vector to reconstructed graph.

        Args:
            z: Latent vector
            data: PyTorch Geometric Data object

        Returns:
            Reconstructed node features
        """
        return self.decoder.decode(z, data)

    def sample(self, num_samples=1, device='cpu'):
        """
        Generate samples from the latent space.

        Args:
            num_samples: Number of samples to generate
            device: Device to generate samples on

        Returns:
            z: Sampled latent vectors
        """
        z = torch.randn(num_samples, self.latent_dim).to(device)
        return z