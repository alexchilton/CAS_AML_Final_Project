import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv, TransformerConv, global_mean_pool


class GATDecoder(nn.Module):
    """
    - Positional embeddings (sequence order)
    - Secondary structure bias (optional)
    """

    def __init__(self, latent_dim, hidden_dim=64, output_dim=None,
                 num_heads=4, dropout=0.1, max_seq_len=3000, ss_dim=8):
        super().__init__()
        self.latent_dim = latent_dim
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim

        # Positional embeddings (learned)
        self.pos_emb = nn.Embedding(max_seq_len, latent_dim)

        # Secondary structure projection (optional)
        self.ss_proj = nn.Linear(ss_dim, latent_dim)

        # Node feature decoder
        self.dec_node_features = nn.Sequential(
            nn.Linear(latent_dim * 2, hidden_dim * 2),  # *2 for pos emb
            nn.ReLU(),
            nn.LayerNorm(hidden_dim * 2),
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.LayerNorm(hidden_dim)
        )

        # GAT layers
        self.gat1 = GATConv(
            hidden_dim, hidden_dim // num_heads,
            heads=num_heads, dropout=dropout
        )
        self.gat2 = GATConv(
            hidden_dim, hidden_dim // num_heads,
            heads=num_heads, dropout=dropout
        )

        # Output projection
        if output_dim is not None:
            self.output_proj = nn.Linear(hidden_dim, output_dim)
        else:
            self.output_proj = None

        self.dropout = nn.Dropout(dropout)

    def forward(self, z, data):
        """Identical signature to original, enhanced with positions/SS"""
        edge_index = data.edge_index
        batch = getattr(data, 'batch', None)
        num_nodes = data.num_nodes

        # Positional embeddings (auto-generated)
        positions = torch.arange(num_nodes, device=z.device)
        pos_emb = self.pos_emb(positions)

        # Combine graph-level z + positional embeddings
        if batch is None:
            graph_emb = z.expand(num_nodes, -1)  # Single graph
        else:
            graph_emb = z[batch]  # Batched graphs
        node_emb = torch.cat([graph_emb, pos_emb], dim=-1)

        # Optional secondary structure bias
        if hasattr(data, 'ss') and data.ss is not None:
            ss_emb = self.ss_proj(data.ss.float())
            node_emb = node_emb[:, :self.latent_dim] + ss_emb

        # Original node feature decoding
        h = self.dec_node_features(node_emb)

        # Original GAT layers
        h = self.gat1(h, edge_index)
        h = F.elu(h)
        h = self.dropout(h)

        h = self.gat2(h, edge_index)
        h = F.elu(h)

        # Dynamic output projection (original behavior)
        if self.output_proj is None:
            self.output_proj = nn.Linear(self.hidden_dim, data.x.size(1)).to(z.device)

        features = self.output_proj(h)
        return self._apply_activations(features)

    def _apply_activations(self, features):
        """Identical to your original activation splitting logic"""
        if features.size(1) >= 39:
            aa_logits = features[:, :21]
            ss_logits = features[:, 21:28]
            coords = features[:, 28:31]
            b_factor = features[:, 31:32]
            meiler = features[:, 32:39]

            aa_probs = F.softmax(aa_logits, dim=1)
            ss_probs = F.softmax(ss_logits, dim=1)
            coords = torch.tanh(coords)
            b_factor = torch.sigmoid(b_factor)
            meiler = torch.sigmoid(meiler)

            return torch.cat([aa_probs, ss_probs, coords, b_factor, meiler,
                              features[:, 39:]], dim=1)
        else:
            return torch.sigmoid(features)

    def decode(self, z, data):
        """Maintains original method signature for compatibility"""
        return self.forward(z, data)