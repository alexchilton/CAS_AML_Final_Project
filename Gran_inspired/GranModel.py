
class GraphAttention(nn.Module):
    def __init__(self, in_features, out_features, n_heads=4, dropout=0.1, alpha=0.2):
        super(GraphAttention, self).__init__()
        self.in_features = in_features
        self.out_features = out_features  # out_features_per_head
        self.n_heads = n_heads
        self.dropout = dropout
        self.alpha = alpha

        # Linear transformations
        self.W = nn.Linear(in_features, n_heads * out_features, bias=False)
        self.a = nn.Linear(2 * out_features, 1, bias=False)

        # Initialize parameters
        nn.init.xavier_uniform_(self.W.weight)
        nn.init.xavier_uniform_(self.a.weight)

        self.leakyrelu = nn.LeakyReLU(self.alpha)
        self.dropout_layer = nn.Dropout(self.dropout)

    def forward(self, h, adj):
        # h: [batch_size, num_nodes, in_features] (32, 50, 128)
        # adj: [batch_size, num_nodes, num_nodes] (32, 50, 50)

        batch_size = h.size(0)
        num_nodes = h.size(1)

        # Linear transformation
        Wh = self.W(h)  # [32, 50, n_heads*out_features] (32,50,4*32=128)
        Wh = Wh.view(batch_size, num_nodes, self.n_heads, self.out_features)  # [32,50,4,32]
        Wh = Wh.permute(0, 2, 1, 3)  # [32,4,50,32]

        # Compute attention coefficients
        Wh_repeated_i = Wh.unsqueeze(3).expand(-1, -1, -1, num_nodes, -1)  # [32,4,50,50,32]
        Wh_repeated_j = Wh.unsqueeze(2).expand(-1, -1, num_nodes, -1, -1)  # [32,4,50,50,32]
        concat = torch.cat([Wh_repeated_i, Wh_repeated_j], dim=-1)  # [32,4,50,50,64]

        e = self.leakyrelu(self.a(concat).squeeze(-1))  # [32,4,50,50]

        # Mask attention coefficients
        zero_vec = -9e15 * torch.ones_like(e)
        adj = adj.unsqueeze(1)  # [32,1,50,50]
        attention = torch.where(adj > 0, e, zero_vec)
        attention = F.softmax(attention, dim=-1)  # [32,4,50,50]
        attention = self.dropout_layer(attention)

        # Apply attention
        h_prime = torch.matmul(attention, Wh)  # [32,4,50,32]
        h_prime = h_prime.permute(0, 2, 1, 3).contiguous()  # [32,50,4,32]
        h_prime = h_prime.view(batch_size, num_nodes, -1)  # [32,50,128]

        return h_prime

#split out
class DualOutputGRAN(nn.Module):
    """
    Graph Recurrent Attention Network that generates both adjacency matrices and amino acid sequences.
    """
    def __init__(self, node_features=38, hidden_dim=128, num_layers=2,
                 n_heads=4, dropout=0.1, amino_acid_vocab_size=22):
        super(DualOutputGRAN, self).__init__()
        self.node_features = node_features
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.n_heads = n_heads
        self.amino_acid_vocab_size = amino_acid_vocab_size

        assert hidden_dim % n_heads == 0, "hidden_dim must be divisible by n_heads"
        self.out_features_per_head = hidden_dim // n_heads

        # Node feature embedding
        self.node_embedding = nn.Linear(node_features, hidden_dim)

        # Graph attention layers
        self.gat_layers = nn.ModuleList()
        for _ in range(num_layers):
            self.gat_layers.append(GraphAttention(
                in_features=hidden_dim,
                out_features=self.out_features_per_head,
                n_heads=n_heads,
                dropout=dropout
            ))

        # Sequence generation
        self.rnn_cell = nn.GRUCell(hidden_dim, hidden_dim)
        self.sequence_projection = nn.Linear(hidden_dim, amino_acid_vocab_size)

        # Secondary structure prediction
        self.ss_projection = nn.Linear(hidden_dim, 9)  # 9 SS types

        # Adjacency matrix generation
        self.edge_predictor = self.create_enhanced_edge_predictor()

        self.dropout = nn.Dropout(dropout)

    def create_enhanced_edge_predictor(self):
        """Enhanced edge predictor with more layers"""
        return nn.Sequential(
            nn.Linear(self.hidden_dim * 2, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(self.hidden_dim, self.hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(self.hidden_dim // 2, 1),
            nn.Sigmoid()
        )

    def compute_advanced_loss(self, predictions, targets):
        """Advanced loss computation with structural constraints and SS prediction"""
        # Existing sequence loss
        seq_logits = predictions['sequence_logits']
        target_seq = targets['sequence']
        batch_size, seq_len, vocab_size = seq_logits.size()
        seq_logits_flat = seq_logits.view(-1, vocab_size)
        target_seq_flat = target_seq.view(-1)
        sequence_loss = F.cross_entropy(seq_logits_flat, target_seq_flat)

        # Enhanced adjacency matrix loss with constraints
        pred_adj = predictions['adjacency_matrix']
        target_adj = targets['adjacency_matrix']

        # 1. Basic contact loss (excluding diagonal)
        mask = 1 - torch.eye(pred_adj.size(1), device=pred_adj.device).unsqueeze(0)
        basic_contact_loss = F.binary_cross_entropy(
            pred_adj * mask,
            target_adj * mask,
            reduction='mean'
        )

        # 2. Sequential distance constraint (consecutive residues ~3.8Å)
        sequential_distance_loss = self.compute_sequential_distance_loss(pred_adj)

        # 3. Symmetry loss
        symmetry_loss = self.compute_symmetry_loss(pred_adj)

        # 4. Distance-based contact classification loss
        distance_loss = self.compute_distance_based_loss(pred_adj, target_adj)

        # Combined adjacency loss with weights
        adjacency_loss = (
                1.0 * basic_contact_loss +
                0.3 * sequential_distance_loss +
                0.3 * symmetry_loss +
                0.4 * distance_loss
        )

        # 5. Secondary structure loss
        ss_logits = predictions['ss_logits']
        # Extract SS labels from node features (positions 7+22 = 29 onwards)
        ss_features = targets['node_features'][:, :, 7+self.amino_acid_vocab_size:]
        ss_labels = torch.argmax(ss_features, dim=-1)
        ss_logits_flat = ss_logits.view(-1, 9)
        ss_labels_flat = ss_labels.view(-1)
        ss_loss = F.cross_entropy(ss_logits_flat, ss_labels_flat)

        # Combined loss with weights
        combined_loss = sequence_loss + adjacency_loss + 0.5 * ss_loss

        return {
            'combined_loss': combined_loss,
            'sequence_loss': sequence_loss,
            'adjacency_loss': adjacency_loss,
            'ss_loss': ss_loss,
            'contact_components': {
                'basic_contact': basic_contact_loss,
                'sequential_distance': sequential_distance_loss,
                'symmetry': symmetry_loss,
                'distance_classification': distance_loss
            }
        }

    def compute_sequential_distance_loss(self, pred_adj):
        """Ensure sequential residues have appropriate distance (~3.8Å)"""
        # Get diagonal bands (distances 1 and 2)
        batch_size, N, _ = pred_adj.size()

        # Distance 1 (adjacent residues)
        diag1 = torch.diagonal(pred_adj, offset=1, dim1=1, dim2=2)
        # Expected probability for distance of 3.8Å in contact maps
        target_prob1 = torch.ones_like(diag1) * 0.9
        loss1 = F.mse_loss(diag1, target_prob1)

        # Distance 2 (residues separated by one)
        diag2 = torch.diagonal(pred_adj, offset=2, dim1=1, dim2=2)
        target_prob2 = torch.ones_like(diag2) * 0.7
        loss2 = F.mse_loss(diag2, target_prob2)

        return (loss1 + loss2) / 2

    def compute_symmetry_loss(self, pred_adj):
        """Enforce symmetry in contact maps"""
        symmetric_diff = torch.abs(pred_adj - pred_adj.transpose(1, 2))
        return torch.mean(symmetric_diff)

    def compute_distance_based_loss(self, pred_adj, target_adj):
        """Classify contacts by distance ranges"""
        batch_size, N, _ = pred_adj.size()

        # Create masks for different distance ranges
        # Local (|i-j| <= 5)
        local_mask = torch.zeros_like(pred_adj)
        for d in range(1, 6):
            local_mask += torch.eye(N, device=pred_adj.device).roll(shifts=d, dims=0).unsqueeze(0)
            local_mask += torch.eye(N, device=pred_adj.device).roll(shifts=-d, dims=0).unsqueeze(0)

        # Medium-range (5 < |i-j| <= 20)
        medium_mask = torch.zeros_like(pred_adj)
        for d in range(6, 21):
            medium_mask += torch.eye(N, device=pred_adj.device).roll(shifts=d, dims=0).unsqueeze(0)
            medium_mask += torch.eye(N, device=pred_adj.device).roll(shifts=-d, dims=0).unsqueeze(0)

        # Long-range (|i-j| > 20)
        long_mask = torch.ones_like(pred_adj) - local_mask - medium_mask
        # Remove diagonal
        long_mask = long_mask * (1 - torch.eye(N, device=pred_adj.device).unsqueeze(0))

        # Compute losses for each range
        local_loss = F.binary_cross_entropy(pred_adj * local_mask, target_adj * local_mask, reduction='sum') / (local_mask.sum() + 1e-8)
        medium_loss = F.binary_cross_entropy(pred_adj * medium_mask, target_adj * medium_mask, reduction='sum') / (medium_mask.sum() + 1e-8)
        long_loss = F.binary_cross_entropy(pred_adj * long_mask, target_adj * long_mask, reduction='sum') / (long_mask.sum() + 1e-8)

        # Weight the losses differently
        return 0.3 * local_loss + 0.3 * medium_loss + 0.4 * long_loss

    def _process_graph(self, node_features, adjacency_matrix):
        """Process the input graph with graph attention layers"""
        h = self.node_embedding(node_features)
        for gat_layer in self.gat_layers:
            h_residual = h  # Save input
            h = gat_layer(h, adjacency_matrix)
            h = F.elu(h + h_residual)  # Add input back to output
            h = self.dropout(h)

        return h

    def _generate_adjacency(self, node_embeddings):
        """Generate adjacency matrix from node embeddings"""
        batch_size, num_nodes, _ = node_embeddings.size()

        # Create all pairwise combinations of node embeddings
        node_i = node_embeddings.unsqueeze(2).repeat(1, 1, num_nodes, 1)  # [B, N, N, H]
        node_j = node_embeddings.unsqueeze(1).repeat(1, num_nodes, 1, 1)  # [B, N, N, H]

        # Concatenate node pairs
        node_pairs = torch.cat([node_i, node_j], dim=-1)  # [B, N, N, 2H]

        # Reshape for passing through edge predictor
        flat_pairs = node_pairs.view(-1, 2 * self.hidden_dim)  # [B*N*N, 2H]

        # Predict edges
        edge_scores = self.edge_predictor(flat_pairs).view(batch_size, num_nodes, num_nodes)  # [B, N, N]

        # Ensure symmetry (for undirected graphs)
        edge_scores = (edge_scores + edge_scores.transpose(1, 2)) / 2

        return edge_scores

    def forward(self, node_features, adjacency_matrix, target_sequences=None, target_adjacency=None, max_length=50):
        """
        Forward pass with dual outputs (sequence and adjacency matrix)

        Args:
            node_features: [batch_size, num_nodes, node_feature_dim]
            adjacency_matrix: [batch_size, num_nodes, num_nodes]
            target_sequences: [batch_size, seq_length] (for training)
            target_adjacency: [batch_size, num_nodes, num_nodes] (for training)
            max_length: Maximum sequence length for generation

        Returns:
            Dictionary containing:
                - 'sequence_logits': Predicted sequence logits
                - 'adjacency_matrix': Predicted adjacency matrix
                - 'ss_logits': Predicted secondary structure logits
                - Or generated sequence and adjacency matrix during inference
        """
        batch_size, num_nodes = node_features.size(0), node_features.size(1)

        # Process graph
        node_embeddings = self._process_graph(node_features, adjacency_matrix)  # [B, N, H]

        # Generate adjacency matrix
        predicted_adjacency = self._generate_adjacency(node_embeddings)  # [B, N, N]

        # Predict secondary structure
        ss_logits = self.ss_projection(node_embeddings)  # [B, N, 9]

        # Initialize RNN for sequence generation
        graph_embedding = torch.mean(node_embeddings, dim=1)  # [B, H]
        h_t = graph_embedding  # Initial hidden state

        if target_sequences is not None:
            # Training mode
            seq_length = target_sequences.size(1)
            sequence_logits = torch.zeros(batch_size, seq_length, self.amino_acid_vocab_size,
                                          device=node_features.device)

            x_t = torch.zeros(batch_size, self.hidden_dim, device=node_features.device)

            for t in range(seq_length):
                h_t = self.rnn_cell(x_t, h_t)  # [B, H]
                sequence_logits[:, t, :] = self.sequence_projection(h_t)  # [B, vocab_size]

                if t < seq_length - 1:
                    # Create a full feature vector for the amino acid
                    amino_acid_feature = torch.zeros(batch_size, self.node_features, device=node_features.device)

                    # Set the amino acid one-hot encoding (positions 7-28)
                    aa_onehot = F.one_hot(target_sequences[:, t], num_classes=self.amino_acid_vocab_size).float()
                    amino_acid_feature[:, 7:7+self.amino_acid_vocab_size] = aa_onehot

                    # Pass through node embedding
                    x_t = self.node_embedding(amino_acid_feature)  # [B, H]

            return {
                'sequence_logits': sequence_logits,
                'adjacency_matrix': predicted_adjacency,
                'ss_logits': ss_logits
            }

        else:
            # Generation mode
            generated_sequences = torch.zeros(batch_size, max_length,
                                              dtype=torch.long,
                                              device=node_features.device)
            x_t = torch.zeros(batch_size, self.hidden_dim,
                              device=node_features.device)

            # Get the full SS predictions upfront
            predicted_ss_full = torch.argmax(ss_logits, dim=-1)  # [B, N]

            for t in range(max_length):
                h_t = self.rnn_cell(x_t, h_t)
                output = self.sequence_projection(h_t)
                prob = F.softmax(output, dim=-1)
                next_aa = torch.multinomial(prob, 1).squeeze(-1)
                generated_sequences[:, t] = next_aa

                # Create a full feature vector for the amino acid
                amino_acid_feature = torch.zeros(batch_size, self.node_features, device=node_features.device)

                # Set the amino acid one-hot encoding (positions 7-28)
                aa_onehot = F.one_hot(next_aa, num_classes=self.amino_acid_vocab_size).float()
                amino_acid_feature[:, 7:7+self.amino_acid_vocab_size] = aa_onehot

                # Add SS prediction if within bounds - simple feedback approach
                if t < predicted_ss_full.size(1):
                    # Get SS prediction for current position
                    ss_onehot = F.one_hot(predicted_ss_full[:, t], num_classes=9).float()
                    # Put SS in positions 29-37
                    amino_acid_feature[:, 7+self.amino_acid_vocab_size:] = ss_onehot

                # Pass through node embedding with both AA and SS info
                x_t = self.node_embedding(amino_acid_feature)  # [B, H]

            return {
                'generated_sequence': generated_sequences,
                'adjacency_matrix': predicted_adjacency,
                'predicted_ss': predicted_ss_full  # [B, N]
            }

