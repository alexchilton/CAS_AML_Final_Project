import pickle
import json

def load_chain_data(chain_output_prefix):
    """
    Load chain data from incremental files.
    
    Parameters:
    -----------
    chain_output_prefix : str
        Prefix used for chain files
        
    Returns:
    --------
    dict
        Combined chain data
    """
    chain_data = {}
    
    try:
        # Try to load index file
        with open(f"{chain_output_prefix}_index.json", 'r') as f:
            index = json.load(f)
        
        # Load each file
        if 'backbone' in index:
            with open(index['backbone'], 'rb') as f:
                chain_data['backbone_atoms'] = pickle.load(f)
        
        if 'secondary_structure' in index:
            with open(index['secondary_structure'], 'rb') as f:
                chain_data['secondary_structure'] = pickle.load(f)
        
        if 'edge_pairs' in index:
            with open(index['edge_pairs'], 'rb') as f:
                chain_data['edge_indices'] = pickle.load(f)
        
        if 'bond_lengths' in index:
            with open(index['bond_lengths'], 'rb') as f:
                chain_data['bond_lengths'] = pickle.load(f)
        
        if 'angles' in index:
            with open(index['angles'], 'rb') as f:
                chain_data['angles'] = pickle.load(f)
        
        if 'torsions' in index:
            with open(index['torsions'], 'rb') as f:
                chain_data['torsions'] = pickle.load(f)
        
        if 'charges' in index:
            with open(index['charges'], 'rb') as f:
                chain_data['charges'] = pickle.load(f)
        
        if 'hydrophobic' in index:
            with open(index['hydrophobic'], 'rb') as f:
                chain_data['hydrophobic_residues'] = pickle.load(f)
        
        return chain_data
        
    except Exception as e:
        print(f"  Error loading chain data: {str(e)}")
        return {}

def visualize_protein_graph(graph_path, output_path=None, show=True, color_by='ss', node_size=100, figsize=(12, 10)):
    """
    Visualize a protein graph from a saved NetworkX pickle file.
    
    Parameters:
    -----------
    graph_path : str
        Path to the saved NetworkX graph pickle file
    output_path : str, optional
        Path to save the visualization image
    show : bool
        Whether to display the graph
    color_by : str
        Property to color nodes by ('ss' for secondary structure,
        'chain', 'residue_name', or any other node attribute)
    node_size : int
        Size of nodes in the visualization
    figsize : tuple
        Figure size (width, height) in inches
    
    Returns:
    --------
    matplotlib figure
    """
    import pickle
    import networkx as nx
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Load the graph
    with open(graph_path, 'rb') as f:
        G = pickle.load(f)
    
    # Create figure
    plt.figure(figsize=figsize)
    
    # Create position layout
    # Use 3D positions if available, otherwise spring layout
    pos = {}
    for node, data in G.nodes(data=True):
        if 'coords' in data:
            # Use first two coordinates for 2D visualization
            pos[node] = data['coords'][:2]
        elif 'CA' in data:
            # If explicit CA coordinates exist
            pos[node] = data['CA'][:2]
    
    # If we couldn't get positions from the graph, use spring layout
    if not pos:
        pos = nx.spring_layout(G, seed=42)
    
    # Color mapping based on the selected property
    if color_by == 'ss':
        # Secondary structure coloring
        ss_types = {'H': 'red', 'E': 'yellow', 'B': 'green', 
                    'G': 'orange', 'I': 'purple', 'T': 'cyan', 
                    'S': 'blue', '?': 'gray', ' ': 'gray'}
        
        node_colors = []
        for node, data in G.nodes(data=True):
            ss = data.get('ss', '?')
            node_colors.append(ss_types.get(ss, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for secondary structure
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=ss) 
                          for ss, color in ss_types.items() if any(data.get('ss', '?') == ss for _, data in G.nodes(data=True))]
        plt.legend(handles=legend_elements, title="Secondary Structure")
        
    elif color_by == 'chain':
        # Color by chain ID
        chains = set()
        for node in G.nodes():
            if ':' in node:
                chains.add(node.split(':')[0])
            else:
                chains.add(G.nodes[node].get('chain_id', 'Unknown'))
        
        # Create color map
        cmap = plt.cm.get_cmap('tab10', len(chains))
        chain_colors = {chain: cmap(i) for i, chain in enumerate(sorted(chains))}
        
        node_colors = []
        for node in G.nodes():
            if ':' in node:
                chain = node.split(':')[0]
            else:
                chain = G.nodes[node].get('chain_id', 'Unknown')
            node_colors.append(chain_colors.get(chain, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for chains
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=chain) 
                          for chain, color in chain_colors.items()]
        plt.legend(handles=legend_elements, title="Chains")
        
    elif color_by == 'residue_name':
        # Color by amino acid type
        aa_types = set()
        for node, data in G.nodes(data=True):
            if ':' in node:
                aa_types.add(node.split(':')[1])
            else:
                aa_types.add(data.get('residue_name', 'Unknown'))
        
        # Create color map
        cmap = plt.cm.get_cmap('tab20', len(aa_types))
        aa_colors = {aa: cmap(i) for i, aa in enumerate(sorted(aa_types))}
        
        node_colors = []
        for node in G.nodes():
            if ':' in node:
                aa = node.split(':')[1]
            else:
                aa = G.nodes[node].get('residue_name', 'Unknown')
            node_colors.append(aa_colors.get(aa, 'gray'))
            
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
        
        # Create legend for amino acids
        from matplotlib.patches import Patch
        # Only show most common AAs in legend to avoid overcrowding
        aa_counts = {}
        for node in G.nodes():
            if ':' in node:
                aa = node.split(':')[1]
            else:
                aa = G.nodes[node].get('residue_name', 'Unknown')
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        # Get top amino acids
        top_aas = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        legend_elements = [Patch(facecolor=aa_colors[aa], label=f"{aa} ({count})") 
                          for aa, count in top_aas]
        plt.legend(handles=legend_elements, title="Top Amino Acids")
    
    else:
        # Try to color by any other attribute
        try:
            values = [G.nodes[node].get(color_by, 0) for node in G.nodes()]
            
            if all(isinstance(v, (int, float)) for v in values if v is not None):
                # Numeric values - use colormap
                cmap = plt.cm.viridis
                vmin = min(v for v in values if v is not None)
                vmax = max(v for v in values if v is not None)
                
                # Replace None with vmin
                values = [vmin if v is None else v for v in values]
                
                nx.draw_networkx_nodes(G, pos, node_color=values, 
                                      cmap=cmap, vmin=vmin, vmax=vmax, node_size=node_size)
                plt.colorbar(label=color_by)
            else:
                # Categorical values
                unique_values = set(v for v in values if v is not None)
                cmap = plt.cm.get_cmap('tab20', len(unique_values))
                value_colors = {val: cmap(i) for i, val in enumerate(sorted(unique_values))}
                
                node_colors = [value_colors.get(G.nodes[node].get(color_by, None), 'gray') 
                              for node in G.nodes()]
                nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_size)
                
                # Create legend
                from matplotlib.patches import Patch
                legend_elements = [Patch(facecolor=color, label=str(val)) 
                                  for val, color in value_colors.items()]
                plt.legend(handles=legend_elements, title=color_by)
        except:
            # Fallback to default coloring
            nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=node_size)
    
    # Draw edges with alpha for clarity
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=0.5)
    
    # Optional: Draw node labels for smaller graphs (uncomment if needed)
    # if G.number_of_nodes() < 100:  # Only show labels for smaller graphs
    #     nx.draw_networkx_labels(G, pos, font_size=8)
    
    # Add graph info
    plt.title(f"Protein Graph: {graph_path.split('/')[-1]}")
    plt.text(0.02, 0.02, f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}", 
             transform=plt.gca().transAxes, fontsize=10)
    
    plt.axis('off')
    
    # Save if requested
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved visualization to {output_path}")
    
    # Show if requested
    if show:
        plt.show()
    else:
        plt.close()
    
    return plt.gcf()

def visualize_all_graphs(output_dir, visualization_dir=None):
    """
    Visualize all graphs in the output directory.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing processed PDB data
    visualization_dir : str, optional
        Directory to save visualizations
    """
    import os
    import glob
    
    if visualization_dir is None:
        visualization_dir = os.path.join(output_dir, "visualizations")
    
    os.makedirs(visualization_dir, exist_ok=True)
    
    # Find all graph files
    graph_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith("_graph.pkl"):
                graph_files.append(os.path.join(root, file))
    
    print(f"Found {len(graph_files)} graph files to visualize")
    
    for graph_file in graph_files:
        pdb_id = os.path.basename(graph_file).replace("_graph.pkl", "")
        print(f"Visualizing {pdb_id}...")
        
        # Create visualizations with different colorings
        for color_scheme in ['ss', 'chain', 'residue_name']:
            output_path = os.path.join(visualization_dir, f"{pdb_id}_{color_scheme}.png")
            
            visualize_protein_graph(
                graph_file,
                output_path=output_path,
                color_by=color_scheme,
                node_size=50,
                show=False  # Don't display, just save
            )
    
    print(f"All visualizations saved to {visualization_dir}")
