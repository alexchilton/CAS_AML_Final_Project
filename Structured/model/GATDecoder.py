import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv
from torch_geometric.data import Data, Batch


class GATDecoder(nn.Module):
    """
    Graph Attention Network (GAT) based decoder for protein graphs.
    Reconstructs a graph from a latent space representation.
    """

    def __init__(self, latent_dim, hidden_dim=64, output_dim=None, num_heads=4, dropout=0.1):
        """
        Initialize the GAT decoder.

        Args:
            latent_dim: Dimension of the latent space
            hidden_dim: Dimension of hidden layers
            output_dim: Dimension of output node features
            num_heads: Number of attention heads for GAT
            dropout: Dropout rate
        """
        super(GATDecoder, self).__init__()

        # Linear layer to transform latent vector to initial node features
        self.fc = nn.Linear(latent_dim, hidden_dim)

        # First GAT layer
        self.gat1 = GATConv(hidden_dim, hidden_dim // num_heads, heads=num_heads, dropout=dropout)

        # Calculate the output dimensions for the second GAT layer
        # We need to ensure the final output has exactly output_dim features
        # If output_dim is not divisible by num_heads, we'll use a different approach
        if output_dim % num_heads == 0:
            # If evenly divisible, we can use output_dim // num_heads for each head
            second_layer_out_dim = output_dim // num_heads
            self.need_final_projection = False
        else:
            # If not evenly divisible, we'll use a standard size and then project
            second_layer_out_dim = hidden_dim // num_heads
            self.need_final_projection = True

        # Second GAT layer
        self.gat2 = GATConv(hidden_dim, second_layer_out_dim, heads=num_heads, dropout=dropout)

        # Final projection layer if needed
        if self.need_final_projection:
            self.final_projection = nn.Linear(second_layer_out_dim * num_heads, output_dim)

        # Dropout layer
        self.dropout = nn.Dropout(dropout)

        # Store output_dim for later reference
        self.output_dim = output_dim

    def forward(self, z, data):
        """
        Forward pass through the decoder.

        Args:
            z: Latent vector
            data: PyTorch Geometric Data object containing the graph structure

        Returns:
            Reconstructed node features
        """
        edge_index = data.edge_index
        batch = data.batch

        # Handle batched graphs
        if batch is None:
            num_nodes = data.x.size(0)
            batch = torch.zeros(num_nodes, dtype=torch.long, device=z.device)
        else:
            # Count nodes per graph in batch
            num_graphs = batch.max().item() + 1
            num_nodes_per_graph = torch.bincount(batch)

        # Transform latent vector to initial node features
        x = self.fc(z)

        # Repeat features for each node in the corresponding graph
        x_list = []
        for i in range(len(z)):
            # For batched input, get number of nodes in this graph
            if batch is not None and i < num_graphs:
                nodes_in_graph = num_nodes_per_graph[i]
                x_repeated = x[i:i+1].repeat(nodes_in_graph, 1)
            else:
                # For single graph
                x_repeated = x[i:i+1].repeat(data.x.size(0), 1)

            x_list.append(x_repeated)

        # Concatenate node features
        x = torch.cat(x_list, dim=0)

        # Apply first GAT layer
        x = self.gat1(x, edge_index)
        x = F.elu(x)
        x = self.dropout(x)

        # Apply second GAT layer
        x = self.gat2(x, edge_index)

        # Apply final projection if needed
        if self.need_final_projection:
            x = self.final_projection(x)

        # Final activation depends on the nature of the features
        # Using sigmoid for normalized features
        x = torch.sigmoid(x)

        return x

    def decode(self, z, data):
        """
        Decode latent vector to reconstructed graph.

        Args:
            z: Latent vector
            data: PyTorch Geometric Data object

        Returns:
            Reconstructed node features
        """
        return self.forward(z, data)