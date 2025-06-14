import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv, global_mean_pool
from torch_geometric.data import Data, Batch


class PositionalEncoding(nn.Module):
    """
    Sinusoidal positional encoding for sequence position information.
    Encodes the sequential ordering of amino acids in protein chains.
    """
    def __init__(self, d_model, max_len=5000):
        """
        Initialize the positional encoding.

        Args:
            d_model: Dimensionality of the model
            max_len: Maximum sequence length to pre-compute
        """
        super(PositionalEncoding, self).__init__()

        # Create constant positional encoding matrix
        # Shape: [max_len, d_model]
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, d_model, 2).float() * (-torch.log(torch.tensor(10000.0)) / d_model))

        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)

        # Register buffer (persistent state)
        self.register_buffer('pe', pe)

    def forward(self, seq_indices):
        """
        Get positional encodings for the given sequence indices.

        Args:
            seq_indices: Tensor of sequence indices [batch_size]

        Returns:
            Positional encodings [batch_size, d_model]
        """
        return self.pe[seq_indices]


class SpatialStructureEmbedding(nn.Module):
    """
    Module for encoding relative spatial positioning of residues using coordinate information.
    Captures the 3D structural relationships between amino acids.
    """
    def __init__(self, input_dim, output_dim):
        """
        Initialize the spatial structure embedding.

        Args:
            input_dim: Input dimension of spatial features
            output_dim: Output dimension of embeddings
        """
        super(SpatialStructureEmbedding, self).__init__()

        # MLP for processing spatial features
        self.mlp = nn.Sequential(
            nn.Linear(input_dim, output_dim * 2),
            nn.ReLU(),
            nn.Linear(output_dim * 2, output_dim)
        )

    def forward(self, coords):
        """
        Process coordinates to extract spatial structure information.

        Args:
            coords: Tensor of 3D coordinates [num_nodes, 3]

        Returns:
            Spatial embeddings [num_nodes, output_dim]
        """
        return self.mlp(coords)


class EnhancedGATEncoder(nn.Module):
    """
    Enhanced Graph Attention Network (GAT) based encoder for protein graphs.
    Incorporates sequence positional embeddings and spatial structure information.
    """

    def __init__(self, input_dim, hidden_dim=64, latent_dim=32, num_heads=4, dropout=0.1,
                 pos_embedding_dim=16, spatial_embedding_dim=16):
        """
        Initialize the enhanced GAT encoder.

        Args:
            input_dim: Dimension of input node features
            hidden_dim: Dimension of hidden layers
            latent_dim: Dimension of the latent space
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
            pos_embedding_dim: Dimension of positional embeddings
            spatial_embedding_dim: Dimension of spatial structure embeddings
        """
        super(EnhancedGATEncoder, self).__init__()

        # Feature dimensions
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim

        # Store indices for sequence position and coordinates in the feature vector
        # Assuming standard feature ordering from ProteinFeatureProcessor:
        # amino acid one-hot (21) + secondary structure one-hot (7) + coordinates (3) + b_factor (1) + meiler features (7)
        self.aa_idx_start = 0
        self.aa_idx_end = 21
        self.coord_idx_start = 21 + 7
        self.coord_idx_end = 21 + 7 + 3

        # Positional encoding for sequence
        self.pos_encoding = PositionalEncoding(d_model=pos_embedding_dim)

        # Spatial structure embedding
        self.spatial_embedding = SpatialStructureEmbedding(
            input_dim=3,  # 3D coordinates
            output_dim=spatial_embedding_dim
        )

        # Feature projection to combine original features with embeddings
        combined_input_dim = input_dim + pos_embedding_dim + spatial_embedding_dim
        self.feature_projection = nn.Linear(combined_input_dim, hidden_dim)

        # First GAT layer
        self.gat1 = GATConv(hidden_dim, hidden_dim // num_heads, heads=num_heads, dropout=dropout)

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
            data: PyTorch Geometric Data object containing the graph and sequence indices

        Returns:
            mu: Mean of the latent distribution
            logvar: Log variance of the latent distribution
        """
        x, edge_index, batch = data.x, data.edge_index, data.batch

        # If batch is None (single graph), create a batch index
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)

        # Extract sequence positions (assume they're stored in the data object)
        # If not available, we can derive them from the edge connectivity
        if hasattr(data, 'seq_pos'):
            seq_positions = data.seq_pos
        else:
            # Default to node indices as a fallback
            seq_positions = torch.arange(x.size(0), device=x.device)

        # Extract coordinates from features
        coords = x[:, self.coord_idx_start:self.coord_idx_end]

        # Get positional encodings
        pos_embeddings = self.pos_encoding(seq_positions)

        # Get spatial structure embeddings
        spatial_embeddings = self.spatial_embedding(coords)

        # Combine features with embeddings
        x_combined = torch.cat([x, pos_embeddings, spatial_embeddings], dim=1)
        x = self.feature_projection(x_combined)

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