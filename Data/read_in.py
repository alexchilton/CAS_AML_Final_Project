import os
import pickle
import json
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def load_protein_data(protein_dir, protein_id=None, chain_id=None):
    """
    Load processed protein data from the output of preprocess_mem_eff.py

    Parameters:
    -----------
    protein_dir : str
        Directory containing the processed data files
    protein_id : str, optional
        Specific protein ID to load (e.g., '1BZQ_nanobody_A')
        If None, will load all proteins in the directory
    chain_id : str, optional
        Specific chain to load (e.g., 'A')

    Returns:
    --------
    dict
        Dictionary of loaded protein data
    """
    result = {}

    # If protein_id is specified, only look for that one
    if protein_id:
        protein_ids = [protein_id]
    else:
        # Get all protein IDs from directories
        protein_ids = []
        for item in os.listdir(protein_dir):
            item_path = os.path.join(protein_dir, item)
            if os.path.isdir(item_path):
                protein_ids.append(item)
            elif item.endswith('_graph.pkl'):
                # Handle files at root level
                protein_ids.append(item.replace('_graph.pkl', ''))

        # Remove duplicates
        protein_ids = list(set(protein_ids))

    for pid in protein_ids:
        protein_path = os.path.join(protein_dir, pid)

        # Initialize data structure
        result[pid] = {}

        # Check if it's a directory structure or flat files
        if os.path.isdir(protein_path):
            # Directory structure

            # Load the main graph
            graph_file = os.path.join(protein_path, f"{pid}_graph.pkl")
            if os.path.exists(graph_file):
                with open(graph_file, 'rb') as f:
                    result[pid]['graph'] = pickle.load(f)
                print(f"Loaded graph for {pid} with {result[pid]['graph'].number_of_nodes()} nodes")

            # Load the main data file
            data_file = os.path.join(protein_path, f"{pid}_data.pkl")
            if os.path.exists(data_file):
                with open(data_file, 'rb') as f:
                    result[pid]['data'] = pickle.load(f)
                print(f"Loaded main data for {pid}")
        else:
            # Flat file structure - files directly in protein_dir

            # Load the main graph
            graph_file = os.path.join(protein_dir, f"{pid}_graph.pkl")
            if os.path.exists(graph_file):
                with open(graph_file, 'rb') as f:
                    result[pid]['graph'] = pickle.load(f)
                print(f"Loaded graph for {pid} with {result[pid]['graph'].number_of_nodes()} nodes")

            # Load the main data file
            data_file = os.path.join(protein_dir, f"{pid}_data.pkl")
            if os.path.exists(data_file):
                with open(data_file, 'rb') as f:
                    result[pid]['data'] = pickle.load(f)
                print(f"Loaded main data for {pid}")

    return result

def visualize_protein_graph(G, output_path=None, color_by='ss'):
    """
    Visualize a protein graph using matplotlib

    Parameters:
    -----------
    G : networkx.Graph
        NetworkX graph to visualize
    output_path : str, optional
        Path to save the visualization. If None, will display instead.
    color_by : str
        Property to color nodes by ('ss', 'chain_id', 'residue_name')
    """
    plt.figure(figsize=(12, 10))

    # Get node positions from coords attribute
    pos = {}
    for node, data in G.nodes(data=True):
        if 'coords' in data:
            # Use only x and y coordinates for 2D visualization
            pos[node] = data['coords'][:2]
        elif 'CA' in data:
            # Use CA atom coordinates if available
            pos[node] = data['CA'][:2]

    # If positions aren't available, use spring layout
    if not pos:
        pos = nx.spring_layout(G, seed=42)

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

        for node, data in G.nodes(data=True):
            ss = data.get('ss', '?')
            node_colors.append(ss_colors.get(ss, 'lightgrey'))

    elif color_by == 'chain_id':
        # Color by chain
        chains = set()
        for _, data in G.nodes(data=True):
            chains.add(data.get('chain_id', 'Unknown'))

        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab10', len(chains))
        chain_colors = {chain: cmap(i) for i, chain in enumerate(sorted(chains))}

        for _, data in G.nodes(data=True):
            chain = data.get('chain_id', 'Unknown')
            node_colors.append(chain_colors.get(chain, 'grey'))

    elif color_by == 'residue_name':
        # Color by amino acid type
        residues = set()
        for _, data in G.nodes(data=True):
            residues.add(data.get('residue_name', 'Unknown'))

        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab20', len(residues))
        residue_colors = {res: cmap(i) for i, res in enumerate(sorted(residues))}

        for _, data in G.nodes(data=True):
            residue = data.get('residue_name', 'Unknown')
            node_colors.append(residue_colors.get(residue, 'grey'))

    else:
        # Default coloring
        node_colors = 'skyblue'

    # Draw the graph
    nx.draw(G, pos, with_labels=False, node_color=node_colors,
            node_size=50, width=0.5, alpha=0.7)

    # Add title
    plt.title("Protein Structure Graph")
    plt.axis('off')

    # Save or display
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
    else:
        plt.show()

def analyze_protein_features(protein_data, protein_id, chain_id):
    """
    Analyze physical features of a protein chain

    Parameters:
    -----------
    protein_data : dict
        Dictionary containing loaded protein data
    protein_id : str
        Identifier for the protein
    chain_id : str
        Chain identifier (e.g., 'A')
    """
    if 'data' not in protein_data[protein_id]:
        print(f"No data found for protein {protein_id}")
        return

    main_data = protein_data[protein_id]['data']

    # Check if chain exists in the main data
    if 'residues_by_chain' in main_data and chain_id not in main_data['residues_by_chain']:
        print(f"Chain {chain_id} not found in protein {protein_id}")
        return

    print(f"\nAnalyzing features for {protein_id}, chain {chain_id}:")

    # Analyze bond lengths
    if 'bond_lengths' in main_data and chain_id in main_data['bond_lengths']:
        bond_lengths = main_data['bond_lengths'][chain_id]

        # Convert the data to a list of values for analysis
        if isinstance(bond_lengths, dict):
            lengths = list(bond_lengths.values())
        else:
            lengths = bond_lengths

        if lengths:
            print(f"\nBond lengths analysis:")
            print(f"  Number of bonds: {len(lengths)}")
            print(f"  Average length: {np.mean(lengths):.3f} Å")
            print(f"  Min length: {min(lengths):.3f} Å")
            print(f"  Max length: {max(lengths):.3f} Å")

    # Analyze angles
    if 'angles' in main_data and chain_id in main_data['angles']:
        angles = main_data['angles'][chain_id]

        # Convert the data to a list of values for analysis
        if isinstance(angles, dict):
            angle_values = list(angles.values())
        else:
            angle_values = angles

        if angle_values:
            print(f"\nBond angles analysis:")
            print(f"  Number of angles: {len(angle_values)}")
            print(f"  Average angle: {np.mean(angle_values):.1f}°")
            print(f"  Min angle: {min(angle_values):.1f}°")
            print(f"  Max angle: {max(angle_values):.1f}°")

    # Analyze torsions
    if 'torsions' in main_data and chain_id in main_data['torsions']:
        torsions = main_data['torsions'][chain_id]

        # Convert the data to a list of values for analysis
        if isinstance(torsions, dict):
            torsion_values = list(torsions.values())
        else:
            torsion_values = torsions

        if torsion_values:
            print(f"\nTorsion angles analysis:")
            print(f"  Number of torsions: {len(torsion_values)}")
            print(f"  Average torsion: {np.mean(torsion_values):.1f}°")
            print(f"  Min torsion: {min(torsion_values):.1f}°")
            print(f"  Max torsion: {max(torsion_values):.1f}°")

    # Analyze secondary structure
    if 'secondary_structure' in main_data and chain_id in main_data['secondary_structure']:
        ss = main_data['secondary_structure'][chain_id]
        ss_counts = {}

        for res_id, ss_type in ss.items():
            ss_counts[ss_type] = ss_counts.get(ss_type, 0) + 1

        if ss_counts:
            print(f"\nSecondary structure analysis:")
            ss_names = {
                'H': 'Alpha helix',
                'B': 'Beta bridge',
                'E': 'Extended strand',
                'G': '3-10 helix',
                'I': 'Pi helix',
                'T': 'Turn',
                'S': 'Bend',
                ' ': 'Coil',
                '?': 'Unknown'
            }
            for ss_type, count in sorted(ss_counts.items()):
                name = ss_names.get(ss_type, ss_type)
                print(f"  {ss_type} ({name}): {count} residues")

    # Analyze hydrophobic residues
    if 'hydrophobic_residues' in main_data and chain_id in main_data['hydrophobic_residues']:
        hydrophobic = main_data['hydrophobic_residues'][chain_id]

        if hydrophobic:
            if isinstance(hydrophobic, set) or isinstance(hydrophobic, list):
                print(f"\nHydrophobic residues: {len(hydrophobic)} residues")
                # Print a few examples
                examples = list(hydrophobic)[:5] if len(hydrophobic) > 5 else list(hydrophobic)
                print(f"  Examples: {examples}")

# Example usage
if __name__ == "__main__":
    # Path to your data
    data_dir = "/Users/alexchilton/Downloads/nanos_networkx"

    # Load the protein and visualize it
    protein_id = "1BZQ_nanobody_A"
    chain_id = "A"

    data = load_protein_data(data_dir, protein_id)

    # Visualize the protein graph with secondary structure coloring
    if protein_id in data and 'graph' in data[protein_id]:
        visualize_protein_graph(data[protein_id]['graph'],
                                output_path=f"{protein_id}_ss.png",
                                color_by='ss')

        # You can also visualize by chain or residue type
        visualize_protein_graph(data[protein_id]['graph'],
                                output_path=f"{protein_id}_chain.png",
                                color_by='chain_id')

        visualize_protein_graph(data[protein_id]['graph'],
                                output_path=f"{protein_id}_residue.png",
                                color_by='residue_name')

    # Analyze features for chain A
    analyze_protein_features(data, protein_id, chain_id)