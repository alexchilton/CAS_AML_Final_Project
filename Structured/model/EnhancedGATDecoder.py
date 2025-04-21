import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv
from torch_geometric.utils import remove_self_loops, add_self_loops
from torch_geometric.data import Data, Batch


class PositionalDecoding(nn.Module):
    """
    Positional decoding module to incorporate sequence information into
    the decoder's node feature generation process.
    """
    def __init__(self, d_model, max_len=5000):
        """
        Initialize the positional decoder.

        Args:
            d_model: Dimensionality of the model
            max_len: Maximum sequence length to pre-compute
        """
        super(PositionalDecoding, self).__init__()

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


class EdgePredictionModule(nn.Module):
    """
    Module for predicting edge connectivity between nodes in the protein graph.
    Evaluates all pairwise connections to determine structure.
    """
    def __init__(self, node_dim, edge_hidden_dim=64, dropout=0.1):
        """
        Initialize the edge prediction module.

        Args:
            node_dim: Dimension of node features
            edge_hidden_dim: Dimension of edge hidden layers
            dropout: Dropout rate
        """
        super(EdgePredictionModule, self).__init__()

        # Edge scoring network
        self.edge_mlp = nn.Sequential(
            nn.Linear(node_dim * 2, edge_hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(edge_hidden_dim, edge_hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(edge_hidden_dim // 2, 1)
        )

        # Threshold for binary edge decisions (can be learned or fixed)
        self.threshold = 0.5

        # Distance-aware edge feature network
        self.edge_feature_network = nn.Sequential(
            nn.Linear(4, edge_hidden_dim),  # 3D distance + edge type
            nn.ReLU(),
            nn.Linear(edge_hidden_dim, 2)   # Edge features
        )

    def forward(self, node_features, num_nodes, temperature=1.0):
        """
        Predict edges based on node features.

        Args:
            node_features: Node feature matrix [num_nodes, node_dim]
            num_nodes: Number of nodes in the graph
            temperature: Temperature for Gumbel-Softmax sampling

        Returns:
            predicted_edge_index: Predicted edge indices [2, num_edges]
            edge_scores: Edge probability scores [num_edges]
            edge_attr: Edge attributes/features [num_edges, edge_attr_dim]
        """
        # Create all possible node pairs (fully connected)
        device = node_features.device

        # Get indices for all unique node pairs
        rows, cols = [], []
        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:  # Exclude self-loops
                    rows.append(i)
                    cols.append(j)

        # Skip if no edges are possible
        if len(rows) == 0:
            return (torch.zeros((2, 0), dtype=torch.long, device=device),
                    torch.zeros((0), device=device),
                    torch.zeros((0, 2), device=device))

        node_pairs_row = torch.tensor(rows, device=device)
        node_pairs_col = torch.tensor(cols, device=device)

        # Get node features for these pairs
        node_features_i = node_features[node_pairs_row]
        node_features_j = node_features[node_pairs_col]

        # Concatenate node features
        edge_features = torch.cat([node_features_i, node_features_j], dim=1)

        # Score edges
        edge_scores = self.edge_mlp(edge_features).squeeze(-1)
        edge_probs = torch.sigmoid(edge_scores)

        # Compute 3D Euclidean distances - assuming 3D coordinates are in node features
        # Extract coordinates (assuming standard ordering from ProteinFeatureProcessor)
        coord_idx_start = 21 + 7  # After one-hot AAs and SS

        # Check if node features has enough dimensions
        if node_features.size(1) > coord_idx_start + 3:
            coords_i = node_features_i[:, coord_idx_start:coord_idx_start+3]
            coords_j = node_features_j[:, coord_idx_start:coord_idx_start+3]
            distances = torch.sqrt(torch.sum((coords_i - coords_j) ** 2, dim=1) + 1e-6)
        else:
            # Fallback if coordinates are not in the expected position
            distances = torch.ones(edge_probs.size(0), device=device)

        # Edge type features (1 for sequential, 0 otherwise - based on node indices)
        sequential_edges = ((node_pairs_row + 1) == node_pairs_col).float()

        # Combine distance with edge type
        edge_input = torch.stack([distances, sequential_edges,
                                  torch.zeros_like(distances),
                                  torch.zeros_like(distances)], dim=1)[:, :4]  # Ensure we have 4 features

        # Apply the Gumbel-Softmax trick for differentiable edge sampling
        if self.training:
            # Sample edges using Gumbel-Softmax
            # logits: [batch_size, 2] - probability of edge existing vs not existing
            logits = torch.stack([edge_probs, 1 - edge_probs], dim=-1)
            edge_decisions = F.gumbel_softmax(logits, tau=temperature, hard=True)[:, 0]
            edge_mask = edge_decisions > 0.5
        else:
            # In evaluation mode, use threshold
            edge_mask = edge_probs > self.threshold

        # Filter selected edges
        selected_indices = torch.nonzero(edge_mask).squeeze(-1)
        if selected_indices.numel() == 0:
            # Handle case with no edges
            return (torch.zeros((2, 0), dtype=torch.long, device=device),
                    torch.zeros((0), device=device),
                    torch.zeros((0, 2), device=device))

        # Create edge index tensor for selected edges
        predicted_edge_index = torch.stack([
            node_pairs_row[selected_indices],
            node_pairs_col[selected_indices]
        ], dim=0)

        # Generate edge attributes for selected edges
        selected_edge_input = edge_input[selected_indices]
        edge_attr = self.edge_feature_network(selected_edge_input)

        # Return edge indices, probabilities, and attributes
        return predicted_edge_index, edge_probs[selected_indices], edge_attr


class EnhancedGATDecoder(nn.Module):
    """
    Enhanced Graph Attention Network (GAT) based decoder for protein graphs.
    Reconstructs node features and predicts edge connectivity simultaneously.
    """

    def __init__(self, latent_dim, hidden_dim=64, output_dim=None, num_heads=4, dropout=0.1,
                 pos_embedding_dim=16, edge_hidden_dim=64):
        """
        Initialize the enhanced GAT decoder.

        Args:
            latent_dim: Dimension of the latent space
            hidden_dim: Dimension of hidden layers
            output_dim: Dimension of output node features
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
            pos_embedding_dim: Dimension of positional embeddings
            edge_hidden_dim: Dimension of edge hidden layers
        """
        super(EnhancedGATDecoder, self).__init__()

        # Store dimensions
        self.latent_dim = latent_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim

        # Positional decoding for sequence
        self.pos_decoding = PositionalDecoding(d_model=pos_embedding_dim)

        # Linear layer to transform latent vector to initial node features
        self.fc_z = nn.Linear(latent_dim, hidden_dim - pos_embedding_dim)

        # MLP for generating spatial coordinates from latent
        self.coord_generator = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 3 * 16)  # Generates coordinates for 16 residues (can be resized)
        )

        # First GAT layer
        self.gat1 = GATConv(hidden_dim, hidden_dim // num_heads, heads=num_heads, dropout=dropout)

        # Calculate the output dimensions for the second GAT layer
        # Ensure the final output has exactly output_dim features
        if output_dim % num_heads == 0:
            # If evenly divisible, use output_dim // num_heads for each head
            second_layer_out_dim = output_dim // num_heads
            self.need_final_projection = False
        else:
            # If not evenly divisible, use a standard size and then project
            second_layer_out_dim = hidden_dim // num_heads
            self.need_final_projection = True

        # Second GAT layer
        self.gat2 = GATConv(hidden_dim, second_layer_out_dim, heads=num_heads, dropout=dropout)

        # Final projection layer if needed
        if self.need_final_projection:
            self.final_projection = nn.Linear(second_layer_out_dim * num_heads, output_dim)

        # Dropout layer
        self.dropout = nn.Dropout(dropout)

        # Edge prediction module
        self.edge_predictor = EdgePredictionModule(
            node_dim=hidden_dim,
            edge_hidden_dim=edge_hidden_dim,
            dropout=dropout
        )

    def generate_node_features(self, z, num_nodes):
        """
        Generate initial node features and positional info from latent vector.

        Args:
            z: Latent vector [batch_size, latent_dim]
            num_nodes: Number of nodes to generate

        Returns:
            node_features: Initial node features [num_nodes, hidden_dim]
            seq_positions: Sequence positions [num_nodes]
        """
        batch_size = z.size(0)
        device = z.device

        # Generate node features from latent vector
        z_features = self.fc_z(z)  # [batch_size, hidden_dim - pos_embedding_dim]

        # Generate sequence positions (consecutive indices)
        seq_positions = torch.arange(num_nodes, device=device)

        # Get positional embeddings
        pos_embeddings = self.pos_decoding(seq_positions)  # [num_nodes, pos_embedding_dim]

        # The key fix: Only use the first graph's latent features in case of batched inputs
        # This ensures consistent dimensions when combining with positional embeddings
        if batch_size > 1:
            z_features = z_features[0:1]  # Use only the first graph's features

        # Expand z features for each node
        z_features_expanded = z_features.repeat(num_nodes, 1)

        # Combine with positional embeddings
        node_features = torch.cat([z_features_expanded, pos_embeddings], dim=1)

        return node_features, seq_positions

    def forward(self, z, num_nodes=None, temperature=1.0):
        """
        Forward pass through the decoder to generate both node features and edges.

        Args:
            z: Latent vector [batch_size, latent_dim]
            num_nodes: Number of nodes to generate (if None, uses a default)
            temperature: Temperature for Gumbel-Softmax edge sampling

        Returns:
            node_features: Reconstructed node features [num_nodes, output_dim]
            edge_index: Predicted edge indices [2, num_edges]
            edge_attr: Edge attributes [num_edges, edge_attr_dim]
        """
        batch_size = z.size(0)
        device = z.device

        # If num_nodes is not provided, estimate it from the latent vector
        # In practice, this would be a learned component or provided externally
        if num_nodes is None:
            # Default to a reasonable size or predict from z
            # For simplicity, using a fixed size here
            num_nodes = 16

        # Generate initial node features and sequence positions
        node_features, seq_positions = self.generate_node_features(z, num_nodes)

        # Predict edge connectivity
        edge_index, edge_probs, edge_attr = self.edge_predictor(
            node_features, num_nodes, temperature
        )

        # If no edges were predicted, create a minimal spanning tree to ensure connectivity
        if edge_index.size(1) == 0:
            # Create a chain-like structure as a fallback
            rows = torch.arange(num_nodes - 1, device=device)
            cols = torch.arange(1, num_nodes, device=device)
            edge_index = torch.stack([rows, cols], dim=0)

            # Create default edge attributes
            edge_attr = torch.ones(edge_index.size(1), 2, device=device)

            # Default edge probabilities
            edge_probs = torch.ones(edge_index.size(1), device=device)

        # Apply the first GAT layer
        x = self.gat1(node_features, edge_index)
        x = F.elu(x)
        x = self.dropout(x)

        # Apply the second GAT layer
        x = self.gat2(x, edge_index)

        # Apply final projection if needed
        if self.need_final_projection:
            x = self.final_projection(x)

        # Final activation depends on the nature of the features
        # Using sigmoid for normalized features
        node_features = torch.sigmoid(x)

        return node_features, edge_index, edge_attr, edge_probs

    def decode(self, z, num_nodes=None):
        """
        Decode latent vector to reconstructed graph (node features and edges).

        Args:
            z: Latent vector [batch_size, latent_dim]
            num_nodes: Number of nodes to generate

        Returns:
            data: PyTorch Geometric Data object with reconstructed graph
        """
        # Generate node features and edges
        node_features, edge_index, edge_attr, _ = self.forward(z, num_nodes)

        # Create a Data object
        data = Data(
            x=node_features,
            edge_index=edge_index,
            edge_attr=edge_attr
        )

        return data