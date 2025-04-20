import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv, global_mean_pool
from torch_geometric.data import Data, Batch


class GATEncoder(nn.Module):
    """
    Graph Attention Network (GAT) based encoder for protein graphs.
    Encodes the input graph into a latent space representation.
    """

    def __init__(self, input_dim, hidden_dim=64, latent_dim=32, num_heads=4, dropout=0.1):
        """
        Initialize the GAT encoder.

        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            latent_dim: Dimension of the latent space
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
        """
        super(GATEncoder, self).__init__()

        # First GAT layer
        self.gat1 = GATConv(input_dim, hidden_dim // num_heads, heads=num_heads, dropout=dropout)

        # Second GAT layer
        self.gat2 = GATConv(hidden_dim, hidden_dim // num_heads, heads=num_heads, dropout=dropout)

        # Linear layers for mean and log variance
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)

        # Dropout layer
        self.dropout = nn.Dropout(dropout)

    def forward(self, data):
        """
        Forward pass through the encoder.

        Args:
            data: PyTorch Geometric Data object containing the graph

        Returns:
            mu: Mean of the latent distribution
            logvar: Log variance of the latent distribution
        """
        x, edge_index, batch = data.x, data.edge_index, data.batch

        # If batch is None (single graph), create a batch index
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        # Apply first GAT layer
        x = self.gat1(x, edge_index)
        x = F.elu(x)
        x = self.dropout(x)

        # Apply second GAT layer
        x = self.gat2(x, edge_index)
        x = F.elu(x)

        # Global mean pooling per graph in batch
        x = global_mean_pool(x, batch)

        # Get latent distribution parameters
        mu = self.fc_mu(x)
        logvar = self.fc_logvar(x)

        return mu, logvar
    def encode(self, data):
        """
        Encode data without returning distribution parameters.
        Useful for getting the latent representation directly.

        Args:
            data: PyTorch Geometric Data object

        Returns:
            z: Sampled latent vector
        """
        mu, logvar = self.forward(data)
        return self.reparameterize(mu, logvar)

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