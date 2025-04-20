import os
import pickle
import json
import numpy as np
import torch
from torch_geometric.data import Data, Dataset
import glob
import matplotlib.pyplot as plt
import networkx as nx

class ProteinGraphDataset(Dataset):
    def __init__(self, protein_dir, transform=None, pre_transform=None, include_features=None):
        """
        Dataset for protein graphs from the preprocess_mem_eff output

        Parameters:
        -----------
        protein_dir : str
            Directory containing the processed data files
        transform : callable, optional
            Transform to apply to each data instance
        pre_transform : callable, optional
            Transform to apply to each data instance before saving
        include_features : list, optional
            List of features to include, e.g., ['bond_lengths', 'angles', 'torsions']
            If None, includes all available features
        """
        super(ProteinGraphDataset, self).__init__(protein_dir, transform, pre_transform)
        self.protein_dir = protein_dir
        self.include_features = include_features or ['bond_lengths', 'angles', 'torsions', 'secondary_structure']

        # Find all available proteins
        self.protein_files = self._find_protein_files()
        print(f"Found {len(self.protein_files)} protein structures")

    def _find_protein_files(self):
        """Find all protein files in the directory"""
        # Look for both graph files and directories
        protein_ids = set()

        # First look for directories
        for item in os.listdir(self.protein_dir):
            item_path = os.path.join(self.protein_dir, item)
            if os.path.isdir(item_path):
                # Check if this directory contains graph and data files
                graph_path = os.path.join(item_path, f"{item}_graph.pkl")
                data_path = os.path.join(item_path, f"{item}_data.pkl")

                if os.path.exists(graph_path) or os.path.exists(data_path):
                    protein_ids.add(item)

        # Also look for graph files directly in the main directory
        graph_files = glob.glob(os.path.join(self.protein_dir, "*_graph.pkl"))
        for graph_file in graph_files:
            protein_id = os.path.basename(graph_file).replace("_graph.pkl", "")
            protein_ids.add(protein_id)

        # Also check for data files
        data_files = glob.glob(os.path.join(self.protein_dir, "*_data.pkl"))
        for data_file in data_files:
            protein_id = os.path.basename(data_file).replace("_data.pkl", "")
            protein_ids.add(protein_id)

        # Convert to list and sort
        protein_ids = sorted(list(protein_ids))

        if protein_ids:
            print(f"Detected proteins: {protein_ids[:5]}... (showing first {min(5, len(protein_ids))})")
        else:
            print("No protein files found! Here are the files in the directory:")
            print(os.listdir(self.protein_dir)[:10])

            # Try a different approach - look for index.json files
            index_files = glob.glob(os.path.join(self.protein_dir, "*_index.json"))
            if index_files:
                print(f"Found {len(index_files)} index files. Sample: {index_files[:3]}")

                # Extract protein IDs from index files
                for idx_file in index_files:
                    base_name = os.path.basename(idx_file)
                    if "_A_index.json" in base_name:
                        protein_id = base_name.replace("_A_index.json", "")
                        protein_ids.append(protein_id)
                    elif "_index.json" in base_name:
                        parts = base_name.split("_")
                        if len(parts) >= 2:
                            # Remove the last two parts ("_index.json")
                            protein_id = "_".join(parts[:-2])
                            protein_ids.append(protein_id)

                protein_ids = sorted(list(set(protein_ids)))
                print(f"Extracted protein IDs: {protein_ids[:5]}... (showing first {min(5, len(protein_ids))})")

        return protein_ids

    def len(self):
        return len(self.protein_files)

    def get(self, idx):
        """Get a single protein graph as a PyG Data object"""
        protein_id = self.protein_files[idx]

        # Load protein data
        data = self._load_protein_data(protein_id)

        if not data or 'graph' not in data:
            print(f"Failed to load valid data for protein {protein_id}")
            # Return an empty Data object as a placeholder
            return Data(x=torch.tensor([[0.0]]), edge_index=torch.zeros((2, 0), dtype=torch.long))

        # Convert to PyG Data
        pyg_data = self._convert_to_pyg(protein_id, data)

        return pyg_data

    def _load_protein_data(self, protein_id):
        """Load data for a specific protein"""
        result = {}

        # Check different possible locations
        possible_locations = [
            # Option 1: Files in a subdirectory named after the protein
            (os.path.join(self.protein_dir, protein_id, f"{protein_id}_graph.pkl"),
             os.path.join(self.protein_dir, protein_id, f"{protein_id}_data.pkl")),

            # Option 2: Files directly in the main directory
            (os.path.join(self.protein_dir, f"{protein_id}_graph.pkl"),
             os.path.join(self.protein_dir, f"{protein_id}_data.pkl"))
        ]

        # Try each location
        for graph_file, data_file in possible_locations:
            # Try to load graph
            if os.path.exists(graph_file):
                try:
                    with open(graph_file, 'rb') as f:
                        result['graph'] = pickle.load(f)
                    print(f"Loaded graph for {protein_id} with {result['graph'].number_of_nodes()} nodes")
                except Exception as e:
                    print(f"Error loading graph for {protein_id}: {e}")

            # Try to load data
            if os.path.exists(data_file):
                try:
                    with open(data_file, 'rb') as f:
                        result['data'] = pickle.load(f)
                    print(f"Loaded data for {protein_id}")
                except Exception as e:
                    print(f"Error loading data for {protein_id}: {e}")

            # If we have at least a graph, we can proceed
            if 'graph' in result:
                break

        # If we still don't have a graph, look for chain-specific files
        if 'graph' not in result:
            print(f"No graph found in standard locations for {protein_id}. Looking for chain-specific files...")

            # Look for index.json files to identify chains
            index_files = glob.glob(os.path.join(self.protein_dir, f"{protein_id}*_index.json"))
            if index_files:
                print(f"Found {len(index_files)} index files for {protein_id}")

                # Try to extract chain information
                chains = []
                for idx_file in index_files:
                    base_name = os.path.basename(idx_file)
                    parts = base_name.split('_')
                    if len(parts) >= 3:
                        chain_id = parts[-2]  # Assuming format like "protein_id_A_index.json"
                        chains.append(chain_id)

                # Try to reconstruct data from chain files
                if chains:
                    print(f"Found chains: {chains}")
                    # For now, just use the first chain to get the graph
                    chain_id = chains[0]

                    # Try to find a backbone file which might contain positions
                    backbone_file = os.path.join(self.protein_dir, f"{protein_id}_{chain_id}_backbone.pkl")
                    if not os.path.exists(backbone_file):
                        backbone_file = os.path.join(self.protein_dir, f"{protein_id}_{chain_id}_{chain_id}_backbone.pkl")

                    if os.path.exists(backbone_file):
                        try:
                            with open(backbone_file, 'rb') as f:
                                backbone_data = pickle.load(f)
                                print(f"Loaded backbone data for {protein_id} chain {chain_id}")

                                # Create a simple graph from backbone data
                                G = nx.Graph()

                                for res_id, atoms in backbone_data.items():
                                    node_id = f"{chain_id}:{res_id}"
                                    G.add_node(node_id,
                                               chain_id=chain_id,
                                               residue_number=res_id,
                                               residue_name="UNK", # We don't have this info
                                               has_backbone=True)

                                    # Add backbone atom coordinates
                                    for atom_name, coords in atoms.items():
                                        G.nodes[node_id][atom_name] = coords
                                        if atom_name == 'CA':
                                            G.nodes[node_id]['coords'] = coords

                                # Add edges between sequential residues
                                res_ids = sorted(backbone_data.keys())
                                for i in range(len(res_ids)-1):
                                    if res_ids[i+1] - res_ids[i] == 1:  # Sequential
                                        G.add_edge(
                                            f"{chain_id}:{res_ids[i]}",
                                            f"{chain_id}:{res_ids[i+1]}",
                                            kind={'peptide_bond'})

                                result['graph'] = G
                                print(f"Created graph with {G.number_of_nodes()} nodes from backbone data")
                        except Exception as e:
                            print(f"Error loading backbone data for {protein_id}: {e}")

        return result

    def _convert_to_pyg(self, protein_id, data):
        """Convert NetworkX graph and associated data to PyG Data object"""
        if 'graph' not in data:
            raise ValueError(f"No graph found for protein {protein_id}")

        nx_graph = data['graph']

        # Get node features
        node_features, node_mapping = self._extract_node_features(nx_graph)

        # Get edge indices
        edge_index, edge_attr = self._extract_edge_features(nx_graph, node_mapping)

        # Extract additional features from data
        additional_features = {}
        if 'data' in data:
            for feature in self.include_features:
                if feature in data['data']:
                    additional_features[feature] = data['data'][feature]

        # Create PyG Data object
        pyg_data = Data(
            x=node_features,
            edge_index=edge_index,
            edge_attr=edge_attr,
            protein_id=protein_id,
            num_nodes=len(node_mapping)
        )

        # Add secondary structure as a node label if available
        if 'data' in data and 'secondary_structure' in data['data']:
            ss_features = self._extract_secondary_structure(nx_graph, node_mapping, data['data']['secondary_structure'])
            if ss_features is not None:
                pyg_data.ss = ss_features

        # Add additional features
        for key, value in additional_features.items():
            if key not in ['secondary_structure']:  # Skip what we've already added
                setattr(pyg_data, key, value)

        return pyg_data

    def _extract_node_features(self, nx_graph):
        """Extract node features from NetworkX graph"""
        # Map node IDs to indices
        nodes = list(nx_graph.nodes())
        node_mapping = {node: i for i, node in enumerate(nodes)}

        # Features we want to extract
        features_to_extract = [
            'chain_id',         # One-hot encoded
            'residue_name',     # One-hot encoded
            'has_backbone',     # Boolean (1/0)
            'coords'            # Coordinates (if available)
        ]

        # Define one-hot encodings
        chain_ids = sorted(set(nx_graph.nodes[node].get('chain_id', 'UNK') for node in nodes))
        chain_id_mapping = {chain: i for i, chain in enumerate(chain_ids)}

        residue_names = sorted(set(nx_graph.nodes[node].get('residue_name', 'UNK') for node in nodes))
        residue_mapping = {res: i for i, res in enumerate(residue_names)}

        # Prepare feature matrix
        feature_dim = len(chain_ids) + len(residue_names) + 1  # chains + residues + has_backbone

        # Add 3 more dimensions for coordinates if they exist
        has_coords = any('coords' in nx_graph.nodes[node] for node in nodes)
        if has_coords:
            feature_dim += 3

        node_features = torch.zeros((len(nodes), feature_dim))

        # Fill feature matrix
        for node, idx in node_mapping.items():
            feature_idx = 0

            # Chain ID (one-hot)
            chain_id = nx_graph.nodes[node].get('chain_id', 'UNK')
            if chain_id in chain_id_mapping:
                node_features[idx, chain_id_mapping[chain_id]] = 1
            feature_idx += len(chain_ids)

            # Residue name (one-hot)
            residue = nx_graph.nodes[node].get('residue_name', 'UNK')
            if residue in residue_mapping:
                node_features[idx, feature_idx + residue_mapping[residue]] = 1
            feature_idx += len(residue_names)

            # Has backbone (binary)
            has_backbone = int(nx_graph.nodes[node].get('has_backbone', False))
            node_features[idx, feature_idx] = has_backbone
            feature_idx += 1

            # Coordinates (if available)
            if has_coords:
                coords = nx_graph.nodes[node].get('coords', None)
                if coords is not None:
                    if isinstance(coords, np.ndarray):
                        coords = torch.from_numpy(coords[:3].copy())  # First 3 dimensions
                    else:
                        coords = torch.tensor(coords[:3])  # First 3 dimensions
                    node_features[idx, feature_idx:feature_idx+3] = coords

        return node_features, node_mapping

    def _extract_edge_features(self, nx_graph, node_mapping):
        """Extract edge indices and features from NetworkX graph"""
        # Collect edges and attributes
        edge_list = []
        edge_attrs = []

        for u, v, data in nx_graph.edges(data=True):
            u_idx = node_mapping[u]
            v_idx = node_mapping[v]

            # Add edge in both directions
            edge_list.append([u_idx, v_idx])
            edge_list.append([v_idx, u_idx])

            # Edge features
            edge_attr = []

            # Edge type (peptide_bond, contact, etc.)
            edge_type = 0  # Default
            if 'kind' in data:
                kind = data['kind']
                if isinstance(kind, set):
                    if 'peptide_bond' in kind:
                        edge_type = 1
                    elif 'contact' in kind:
                        edge_type = 2
                elif isinstance(kind, str):
                    if kind == 'peptide_bond':
                        edge_type = 1
                    elif kind == 'contact':
                        edge_type = 2

            edge_attr.append(edge_type)

            # Distance (if available)
            if 'distance' in data:
                edge_attr.append(data['distance'])
            else:
                edge_attr.append(0.0)  # Default

            # Add edge attributes (both directions)
            edge_attrs.append(edge_attr)
            edge_attrs.append(edge_attr)

        # Convert to tensors
        if edge_list:
            edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
            edge_attr = torch.tensor(edge_attrs, dtype=torch.float)
        else:
            # Empty graph
            edge_index = torch.zeros((2, 0), dtype=torch.long)
            edge_attr = torch.zeros((0, 2), dtype=torch.float)

        return edge_index, edge_attr

    def _extract_secondary_structure(self, nx_graph, node_mapping, ss_data):
        """Extract secondary structure labels from data"""
        if not ss_data:
            return None

        # SS types
        ss_types = {
            'H': 0,  # Alpha helix
            'B': 1,  # Beta bridge
            'E': 2,  # Extended strand
            'G': 3,  # 3-10 helix
            'I': 4,  # Pi helix
            'T': 5,  # Turn
            'S': 6,  # Bend
            ' ': 7,  # Coil
            '?': 8   # Unknown
        }

        ss_features = torch.zeros(len(node_mapping), dtype=torch.long)

        # Map secondary structure to nodes
        for node, idx in node_mapping.items():
            # Extract chain and residue ID
            node_data = nx_graph.nodes[node]
            chain_id = node_data.get('chain_id', None)

            if chain_id and chain_id in ss_data:
                # Get residue ID
                if ':' in node:
                    # Format like 'A:LYS:1'
                    res_id = int(node.split(':')[-1])
                else:
                    # Get from attributes
                    res_id = node_data.get('residue_number', None)

                if res_id and res_id in ss_data[chain_id]:
                    ss_type = ss_data[chain_id][res_id]
                    ss_features[idx] = ss_types.get(ss_type, 8)  # Default to unknown

        return ss_features

def visualize_protein_graph(graph, protein_id=None, output_path=None, color_by='ss'):
    """Visualize a protein graph using NetworkX and matplotlib"""
    plt.figure(figsize=(12, 10))

    # Get node positions from coords attribute
    pos = {}
    for node, data in graph.nodes(data=True):
        if 'coords' in data:
            # Use only x and y coordinates for 2D visualization
            pos[node] = data['coords'][:2]
        elif 'CA' in data:
            # Use CA atom coordinates if available
            pos[node] = data['CA'][:2]

    # If positions aren't available, use spring layout
    if not pos:
        pos = nx.spring_layout(graph, seed=42)

    # Color nodes based on selected property
    node_colors = []

    if color_by == 'ss':
        # Color by secondary structure
        ss_colors = {
            'H': 'red',      # Alpha helix
            'B': 'blue',     # Beta bridge
            'E': 'yellow',   # Extended strand
            'G': 'orange',   # 3-10 helix
            'I': 'purple',   # Pi helix
            'T': 'cyan',     # Turn
            'S': 'green',    # Bend
            ' ': 'grey',     # Loop/irregular
            '?': 'lightgrey' # Unknown
        }

        for node, data in graph.nodes(data=True):
            ss = data.get('ss', '?')
            node_colors.append(ss_colors.get(ss, 'lightgrey'))

    elif color_by == 'chain_id':
        # Color by chain
        chains = set()
        for _, data in graph.nodes(data=True):
            chains.add(data.get('chain_id', 'Unknown'))

        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab10', len(chains))
        chain_colors = {chain: cmap(i) for i, chain in enumerate(sorted(chains))}

        for _, data in graph.nodes(data=True):
            chain = data.get('chain_id', 'Unknown')
            node_colors.append(chain_colors.get(chain, 'grey'))

    elif color_by == 'residue_name':
        # Color by amino acid type
        residues = set()
        for _, data in graph.nodes(data=True):
            residues.add(data.get('residue_name', 'Unknown'))

        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab20', len(residues))
        residue_colors = {res: cmap(i) for i, res in enumerate(sorted(residues))}

        for _, data in graph.nodes(data=True):
            residue = data.get('residue_name', 'Unknown')
            node_colors.append(residue_colors.get(residue, 'grey'))

    else:
        # Default coloring
        node_colors = 'skyblue'

    # Draw the graph
    nx.draw(graph, pos, with_labels=False, node_color=node_colors,
            node_size=50, width=0.5, alpha=0.7)

    # Add title
    title = f"Protein Structure Graph: {protein_id}" if protein_id else "Protein Structure Graph"
    plt.title(title)
    plt.axis('off')

    # Save or display
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
    else:
        plt.show()

def load_and_convert_all_proteins(protein_dir, output_dir=None, visualize=False):
    """
    Load all proteins and convert them to PyTorch Geometric format

    Parameters:
    -----------
    protein_dir : str
        Directory containing the processed protein data
    output_dir : str, optional
        Directory to save the converted data
    visualize : bool
        Whether to visualize the first few proteins

    Returns:
    --------
    list
        List of PyG Data objects
    """
    # Create dataset
    dataset = ProteinGraphDataset(protein_dir)

    if len(dataset) == 0:
        print("No proteins found in dataset!")
        return []

    # Process and optionally save all proteins
    all_data = []
    for i in range(len(dataset)):
        protein_id = dataset.protein_files[i]
        print(f"Processing protein {i+1}/{len(dataset)}: {protein_id}")

        try:
            # Get PyG data
            pyg_data = dataset[i]
            all_data.append(pyg_data)

            # Visualize if requested
            if visualize and i < 5:  # Visualize first 5 proteins
                # Load original graph for visualization
                protein_data = dataset._load_protein_data(protein_id)
                if 'graph' in protein_data:
                    visualize_protein_graph(
                        protein_data['graph'],
                        protein_id=protein_id,
                        output_path=os.path.join(output_dir, f"{protein_id}_graph.png") if output_dir else None,
                        color_by='ss'
                    )

            # Save to disk if output directory provided
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                torch.save(pyg_data, os.path.join(output_dir, f"{protein_id}.pt"))
        except Exception as e:
            print(f"Error processing protein {protein_id}: {e}")

    print(f"Processed {len(all_data)} proteins")
    return all_data

def create_gnn_dataset_from_directory(input_dir, output_dir=None, visualize=False):
    """
    Main function to create a PyTorch Geometric dataset from preprocessed protein data

    Parameters:
    -----------
    input_dir : str
        Directory containing preprocessed protein data
    output_dir : str, optional
        Directory to save the PyG data and visualizations
    visualize : bool
        Whether to visualize some examples

    Returns:
    --------
    list
        List of PyG Data objects
    """
    # Create output directory if needed
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Load and convert all proteins
    all_data = load_and_convert_all_proteins(input_dir, output_dir, visualize)

    # Save the whole dataset
    if output_dir and all_data:
        torch.save(all_data, os.path.join(output_dir, "protein_dataset.pt"))

    return all_data

# Example usage
if __name__ == "__main__":
    input_dir = "/Users/alexchilton/Downloads/nanos_networkx"
    output_dir = "/Users/alexchilton/Downloads/nanos_pyg"

    # Create PyG dataset
    protein_dataset = create_gnn_dataset_from_directory(input_dir, output_dir, visualize=True)

    # Print some statistics
    print("\nDataset overview:")
    print(f"Number of proteins: {len(protein_dataset)}")

    # Print properties of first protein
    if protein_dataset:
        print("\nFirst protein properties:")
        data = protein_dataset[0]
        print(f"Protein ID: {data.protein_id}")
        print(f"Nodes: {data.num_nodes}")
        print(f"Edges: {data.edge_index.size(1)}")
        print(f"Node features: {data.x.size()}")
        print(f"Edge features: {data.edge_attr.size()}")
        if hasattr(data, 'ss'):
            print(f"Secondary structure labels: {data.ss.size()}")