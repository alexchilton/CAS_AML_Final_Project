
def load_protein_graph_data(parent_folder, max_proteins=None, chunk_length=chunk_length):
    """
    Load protein NetworkX graph data with detailed diagnostics
    """
    random.seed(42)

    # List to store graphs and their associated amino acid sequences
    protein_graphs = []
    protein_sequences = []

    # Diagnostic counters
    loaded_proteins = 0
    failed_proteins = 0

    folders = [name for name in os.listdir(parent_folder)
               if os.path.isdir(os.path.join(parent_folder, name))]

    found_folders = len(folders)
    print(f"Found {found_folders} protein folders")

    # Load only max_proteins if specified
    if max_proteins:
        folders = folders[:max_proteins]
        print(f"Limited to {len(folders)} folders due to max_proteins setting")

    for folder in folders:
        # Look for graph files in the folder
        folder_path = os.path.join(parent_folder, folder)
        graph_files = [f for f in os.listdir(folder_path) if f.endswith('_graph.pkl')]

        if graph_files:
            graph_file = os.path.join(folder_path, graph_files[0])
            try:
                with open(graph_file, 'rb') as f:
                    graph = pkl.load(f)

                    # Extract amino acid sequence
                    aa_seq = []
                    sorted_nodes = sorted(graph.nodes(), key=lambda x:
                    int(graph.nodes[x]['residue_number'])
                    if 'residue_number' in graph.nodes[x] else 0)

                    for node in sorted_nodes:
                        if 'residue_name' in graph.nodes[node]:
                            aa_seq.append(graph.nodes[node]['residue_name'])
                        else:
                            aa_seq.append('X')

                    protein_graphs.append(graph)
                    protein_sequences.append(aa_seq)
                    loaded_proteins += 1
            except Exception as e:
                failed_proteins += 1
                print(f"Error loading {graph_file}: {e}")

    print(f"Successfully loaded {loaded_proteins} protein graphs, {failed_proteins} failed")

    # Create smaller subgraphs
    subgraphs = []
    subsequences = []

    # Diagnostic counters for subsequence generation
    total_possible_subsequences = 0
    actual_generated_subsequences = 0
    proteins_with_no_subsequences = 0

    for i, (graph, sequence) in enumerate(zip(protein_graphs, protein_sequences)):
        # Extract nodes by residue number ranges
        sorted_nodes = sorted(graph.nodes(), key=lambda x:
        int(graph.nodes[x]['residue_number'])
        if 'residue_number' in graph.nodes[x] else 0)

        # Count theoretical maximum for this protein
        seq_length = len(sorted_nodes)
        max_subseqs_this_protein = seq_length  # With wrapping, we can generate as many subsequences as residues
        total_possible_subsequences += max_subseqs_this_protein

        subseqs_this_protein = 0

        # Create overlapping subsequences with wrapping
        for start_idx in range(0, seq_length, step_size):
            # Generate indices for the subsequence with wrapping
            indices = []
            for j in range(chunk_length):
                # Wrap around if we exceed the length
                wrapped_idx = (start_idx + j) % seq_length
                indices.append(wrapped_idx)

            # Get the corresponding nodes
            node_subset = [sorted_nodes[idx] for idx in indices]
            # Get the corresponding sequence
            aa_subset = [sequence[idx] for idx in indices]

            if node_subset:
                try:
                    subgraph = graph.subgraph(node_subset)
                    if len(subgraph) > 0:
                        subgraphs.append(subgraph)
                        subsequences.append(aa_subset)
                        subseqs_this_protein += 1
                        actual_generated_subsequences += 1
                except Exception as e:
                    print(f"Error creating subgraph: {e}")

        if subseqs_this_protein == 0:
            proteins_with_no_subsequences += 1

    print(f"\nSubsequence Generation Statistics:")
    print(f"Total possible subsequences: {total_possible_subsequences}")
    print(f"Actually generated subsequences: {actual_generated_subsequences}")
    print(f"Proteins with no subsequences: {proteins_with_no_subsequences}")

    # If no subgraphs were created, use the full graphs
    if not subgraphs:
        print("Warning: Could not create subgraphs. Using full graphs instead.")
        subgraphs = protein_graphs
        subsequences = protein_sequences

    print(f"Final count: {len(subgraphs)} protein subgraphs for training")

    # Shuffle the data
    if len(subgraphs) > 0:
        combined = list(zip(subgraphs, subsequences))
        random.shuffle(combined)
        subgraphs, subsequences = zip(*combined)

    # Limit size if needed
    if max_proteins and len(subgraphs) > subgraph_limit:
        print(f"Limiting from {len(subgraphs)} to {subgraph_limit} due to max_proteins /subgraph limit")
        subgraphs = subgraphs[:subgraph_limit]
        subsequences = subsequences[:subgraph_limit]

    return protein_graphs, protein_sequences, list(subgraphs), list(subsequences)

def prepare_graph_data_for_training(subgraphs, subsequences, unique_aa):
    """
    Prepare graph data for GRAN model training with enhanced node features

    Args:
        subgraphs: List of NetworkX subgraphs
        subsequences: List of amino acid sequences corresponding to the subgraphs
        unique_aa: Set of unique amino acids to use for one-hot encoding

    Returns:
        aa_sequences_tensor: Tensor of amino acid sequences (targets)
        adjacency_tensors: Tensor of graph adjacency matrices
        node_features_tensor: Tensor of node features
    """
    # Secondary structure mapping
    SS_TO_INDEX = {
        'E': 0,  # Extended (beta sheets)
        '-': 1,  # Coil
        'T': 2,  # Turn
        'S': 3,  # Bend
        'G': 4,  # 3-10 helix
        'H': 5,  # Alpha helix
        'B': 6,  # Bridge
        'I': 7,  # Pi helix
        '?': 8   # Unknown
    }

    # Prepare amino acid sequences (these will be our targets)
    aa_sequences = []
    for seq in subsequences:
        # Convert amino acid sequence to indices
        aa_indices = []
        for aa in seq:
            if aa in unique_aa:
                aa_indices.append(list(unique_aa).index(aa))
            else:
                # Handle unknown amino acids
                aa_indices.append(list(unique_aa).index('X') if 'X' in unique_aa else 0)
        aa_sequences.append(aa_indices)

    # Prepare adjacency matrices and node features
    adjacency_matrices = []
    node_features = []

    for i, graph in enumerate(subgraphs):
        # Get number of nodes
        num_nodes = len(graph)

        # Skip empty graphs
        if num_nodes == 0:
            continue

        # Create adjacency matrix
        adj_matrix = nx.to_numpy_array(graph)
        adjacency_matrices.append(adj_matrix)

        # Extract ordered nodes
        sorted_nodes = sorted(graph.nodes(), key=lambda x:
        int(graph.nodes[x]['residue_number'])
        if 'residue_number' in graph.nodes[x] else 0)

        # Create enhanced features: Meiler (7) + AA one-hot (22) + SS one-hot (9) = 38
        FEATURE_DIM = 7 + len(unique_aa) + 9
        features = np.zeros((num_nodes, FEATURE_DIM))

        for j, node in enumerate(sorted_nodes):
            # Meiler features (first 7 dimensions)
            if 'meiler' in graph.nodes[node]:
                meiler_values = graph.nodes[node]['meiler']
                # Extract meiler dimensions - handle different formats
                if isinstance(meiler_values, dict):
                    # Format: {'dim_1': 1.28, 'dim_2': 0.05, ...}
                    features[j, 0] = meiler_values.get('dim_1', 0.0)
                    features[j, 1] = meiler_values.get('dim_2', 0.0)
                    features[j, 2] = meiler_values.get('dim_3', 0.0)
                    features[j, 3] = meiler_values.get('dim_4', 0.0)
                    features[j, 4] = meiler_values.get('dim_5', 0.0)
                    features[j, 5] = meiler_values.get('dim_6', 0.0)
                    features[j, 6] = meiler_values.get('dim_7', 0.0)
                elif isinstance(meiler_values, list) and len(meiler_values) >= 7:
                    # Format: [1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23]
                    features[j, :7] = meiler_values[:7]

            # Amino acid one-hot (dimensions 7-28)
            if 'residue_name' in graph.nodes[node]:
                aa = graph.nodes[node]['residue_name']
                if aa in unique_aa:
                    aa_idx = list(unique_aa).index(aa)
                    features[j, 7 + aa_idx] = 1.0
                elif 'X' in unique_aa:
                    aa_idx = list(unique_aa).index('X')
                    features[j, 7 + aa_idx] = 1.0
                else:
                    features[j, 7] = 1.0

            # Secondary structure one-hot (dimensions 29-37)
            if 'secondary_structure' in graph.nodes[node]:
                ss = graph.nodes[node]['secondary_structure']
                if ss in SS_TO_INDEX:
                    ss_idx = SS_TO_INDEX[ss]
                    features[j, 7 + len(unique_aa) + ss_idx] = 1.0
                else:
                    # Unknown SS -> use '?' index
                    features[j, 7 + len(unique_aa) + SS_TO_INDEX['?']] = 1.0

        node_features.append(features)

    # Check if we have data
    if not adjacency_matrices or not node_features:
        raise ValueError("No valid graphs found after processing")

    # Convert to tensors
    aa_sequences_tensor = [torch.tensor(seq, dtype=torch.long) for seq in aa_sequences if seq]
    adjacency_tensors = [torch.tensor(adj, dtype=torch.float32) for adj in adjacency_matrices]
    node_features_tensor = [torch.tensor(nf, dtype=torch.float32) for nf in node_features]

    return aa_sequences_tensor, adjacency_tensors, node_features_tensor

def collate_batch(batch):
    aa_seqs, adjacency_matrices, node_feats = zip(*batch)
    # Stack directly without padding since all should be size 50
    return torch.stack(aa_seqs), torch.stack(adjacency_matrices), torch.stack(node_feats)


def create_dataloader(aa_sequences, adjacency_matrices, node_features, batch_size=batch_size, shuffle=True):
    """
    Create a dataloader from the prepared data
    """
    dataset = list(zip(aa_sequences, adjacency_matrices, node_features))
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        collate_fn=collate_batch
    )
    return dataloader

