import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data

from EnhancedGATEncoder import EnhancedGATEncoder
from EnhancedGATDecoder import EnhancedGATDecoder


class EnhancedGATVAE(nn.Module):
    """
    Enhanced Graph Attention Variational Autoencoder for protein graphs.
    Combines the EnhancedGATEncoder and EnhancedGATDecoder to create a model
    that can learn both node features and graph structure.
    """

    def __init__(self, input_dim, hidden_dim=64, latent_dim=32, num_heads=4, dropout=0.1,
                 pos_embedding_dim=16, spatial_embedding_dim=16, edge_hidden_dim=64):
        """
        Initialize the enhanced GAT-based VAE.

        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            latent_dim: Dimension of the latent space
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
            pos_embedding_dim: Dimension of positional embeddings
            spatial_embedding_dim: Dimension of spatial structure embeddings
            edge_hidden_dim: Dimension of edge hidden layers
        """
        super(EnhancedGATVAE, self).__init__()

        # Initialize enhanced encoder
        self.encoder = EnhancedGATEncoder(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            latent_dim=latent_dim,
            num_heads=num_heads,
            dropout=dropout,
            pos_embedding_dim=pos_embedding_dim,
            spatial_embedding_dim=spatial_embedding_dim
        )

        # Initialize enhanced decoder
        self.decoder = EnhancedGATDecoder(
            latent_dim=latent_dim,
            hidden_dim=hidden_dim,
            output_dim=input_dim,
            num_heads=num_heads,
            dropout=dropout,
            pos_embedding_dim=pos_embedding_dim,
            edge_hidden_dim=edge_hidden_dim
        )

        self.latent_dim = latent_dim

    def forward(self, data):
        """
        Forward pass through the VAE. Encodes input and then decodes to reconstruct
        both node features and graph structure.

        Args:
            data: PyTorch Geometric Data object containing the graph

        Returns:
            reconstructed_data: Data object with reconstructed features and structure
            mu: Mean of the latent distribution
            logvar: Log variance of the latent distribution
            edge_probabilities: Probabilities for predicted edges
        """
        # Encode
        mu, logvar = self.encoder(data)

        # Sample from latent distribution
        z = self.encoder.reparameterize(mu, logvar)

        # Decode to get node features and graph structure
        # Use the number of nodes from the input data
        num_nodes = data.x.size(0)
        node_features, edge_index, edge_attr, edge_probs = self.decoder(z, num_nodes)

        # Create reconstructed data object
        reconstructed_data = Data(
            x=node_features,
            edge_index=edge_index,
            edge_attr=edge_attr
        )

        return reconstructed_data, mu, logvar, edge_probs

    def encode(self, data):
        """
        Encode data to latent representation.

        Args:
            data: PyTorch Geometric Data object

        Returns:
            z: Latent vector
        """
        return self.encoder.encode(data)

    def decode(self, z, num_nodes=None):
        """
        Decode latent vector to reconstructed graph.

        Args:
            z: Latent vector
            num_nodes: Number of nodes to generate (if None, model estimates it)

        Returns:
            Reconstructed graph as a Data object
        """
        return self.decoder.decode(z, num_nodes)

    def sample(self, num_samples=1, num_nodes=16, device='cpu'):
        """
        Generate new protein structures by sampling from the latent space.

        Args:
            num_samples: Number of samples to generate
            num_nodes: Number of nodes for each generated protein
            device: Device to generate samples on

        Returns:
            List of sampled graphs as Data objects
        """
        z = torch.randn(num_samples, self.latent_dim, device=device)

        samples = []
        for i in range(num_samples):
            # Generate a single sample
            sample = self.decode(z[i:i+1], num_nodes)
            samples.append(sample)

        return samples