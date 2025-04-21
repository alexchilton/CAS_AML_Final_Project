import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data

from EnhancedGATEncoder import EnhancedGATEncoder
from FixedEnhancedGATDecoder import FixedEnhancedGATDecoder

class FixedEnhancedGATVAE(nn.Module):
    """
    Fixed version of the Enhanced Graph Attention Variational Autoencoder for protein graphs.
    Uses the improved decoder that generates diverse node features.
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
        super(FixedEnhancedGATVAE, self).__init__()

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

        # Initialize fixed enhanced decoder
        self.decoder = FixedEnhancedGATDecoder(
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
        z = self.reparameterize(mu, logvar)

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
        mu, logvar = self.encoder(data)
        return self.reparameterize(mu, logvar)

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

    def reparameterize(self, mu, logvar):
        """
        Reparameterization trick to sample from N(mu, var) from N(0,1).

        Args:
            mu: Mean of the latent Gaussian
            logvar: Log variance of the latent Gaussian

        Returns:
            Sampled latent vector
        """
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def sample(self, num_samples=1, num_nodes=16, device='cpu', dataset_handler=None):
        """
        Generate new protein structures by sampling from the latent space.

        FIXED: Now ensures different seeds for each sample.

        Args:
            num_samples: Number of samples to generate
            num_nodes: Number of nodes for each generated protein
            device: Device to generate samples on
            dataset_handler: Optional dataset handler for denormalizing data

        Returns:
            List of sampled graphs as Data objects
        """
        samples = []

        for i in range(num_samples):
            # Use a different seed for each sample to ensure diversity
            torch.manual_seed(42 + i)

            # Create random latent vector
            z = torch.randn(1, self.latent_dim, device=device)

            # Generate a single sample
            sample = self.decode(z, num_nodes)

            # Denormalize if dataset_handler is provided
            if dataset_handler is not None:
                sample = dataset_handler.denormalize_data(sample)

                # Optionally convert continuous AA features to discrete one-hot encoding
                if hasattr(sample, 'x') and sample.x.size(1) >= 21:  # Check if we have AA features
                    aa_features = sample.x[:, :21]
                    aa_indices = torch.argmax(aa_features, dim=1)

                    # Create new one-hot encodings
                    new_aa_features = torch.zeros_like(aa_features)
                    for j, idx in enumerate(aa_indices):
                        new_aa_features[j, idx] = 1.0

                    # Replace the original features
                    sample.x[:, :21] = new_aa_features

            samples.append(sample)

        return samples