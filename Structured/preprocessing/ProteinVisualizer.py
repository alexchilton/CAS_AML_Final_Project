import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

class ProteinVisualizer:
    """
    A class for visualizing protein structures and their graph representations.
    
    This class handles:
    - NetworkX graph visualization
    - Graph inspection and property reporting
    - Structure visualization
    """

    def __init__(self, output_dir=None):
        """
        Initialize the protein visualizer.
        
        Parameters:
        -----------
        output_dir : str, optional
            Directory to save visualizations. If None, visualizations are displayed
            but not saved.
        """
        self.output_dir = output_dir
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    def visualize_protein_graph(self, graph, pdb_id=None, output_path=None, figsize=(12, 10)):
        """
        Visualize the protein graph using matplotlib.
        
        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to visualize
        pdb_id : str, optional
            PDB identifier for title
        output_path : str, optional
            Path to save the visualization
        figsize : tuple, optional
            Figure size as (width, height)
            
        Returns:
        --------
        None
            Displays or saves the plot
        """
        if graph.number_of_nodes() == 0:
            print("Cannot visualize empty graph!")
            return

        # Set up figure
        plt.figure(figsize=figsize)

        # Collect node positions and colors
        pos = {}
        colors = []
        sizes = []

        # Determine if 'coords' attribute exists or needs to be derived
        has_coords = False
        for node in list(graph.nodes())[:5]:  # Check first few nodes
            if 'coords' in graph.nodes[node]:
                has_coords = True
                break

        # For each node, get position and color
        for node in graph.nodes():
            # Get chain_id
            chain_id = graph.nodes[node].get('chain_id', 'Unknown')

            # Get position
            if has_coords and 'coords' in graph.nodes[node]:
                coords = graph.nodes[node]['coords']
                if isinstance(coords, np.ndarray) and coords.size >= 2:
                    pos[node] = (coords[0], coords[1])
                elif isinstance(coords, list) and len(coords) >= 2:
                    pos[node] = (coords[0], coords[1])
            elif 'CA' in graph.nodes[node]:
                # Use CA atom position if coords not available
                ca_pos = graph.nodes[node]['CA']
                if isinstance(ca_pos, np.ndarray) and ca_pos.size >= 2:
                    pos[node] = (ca_pos[0], ca_pos[1])
                elif isinstance(ca_pos, list) and len(ca_pos) >= 2:
                    pos[node] = (ca_pos[0], ca_pos[1])

            # Add color based on chain
            if chain_id == 'L':
                colors.append('blue')
            elif chain_id == 'H':
                colors.append('red')
            else:
                colors.append('gray')

            # Node size
            sizes.append(30)  # Default size

        # If we couldn't get positions for nodes, use spring layout
        if len(pos) < graph.number_of_nodes() * 0.5:  # If less than half the nodes have positions
            print("Using spring layout as fallback")
            pos = nx.spring_layout(graph, seed=42)

        # Draw the graph
        nx.draw(
            graph,
            pos=pos,
            node_color=colors,
            node_size=sizes,
            with_labels=False,
            edge_color='gray',
            width=0.5,
            alpha=0.7
        )

        # Add title
        title = f"Protein Graph: {pdb_id}" if pdb_id else "Protein Graph"
        plt.title(title)
        plt.axis('off')

        # Add legend
        legend_handles = []

        if 'blue' in colors:
            legend_handles.append(mlines.Line2D([], [], color='blue', marker='o',
                                                linestyle='None', markersize=10, label='Chain L'))
        if 'red' in colors:
            legend_handles.append(mlines.Line2D([], [], color='red', marker='o',
                                                linestyle='None', markersize=10, label='Chain H'))
        if 'gray' in colors:
            legend_handles.append(mlines.Line2D([], [], color='gray', marker='o',
                                                linestyle='None', markersize=10, label='Other'))

        plt.legend(handles=legend_handles)

        # Save or show
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def visualize_protein_graph_by_ss(self, graph, pdb_id=None, output_path=None, figsize=(12, 10)):
        """
        Visualize the protein graph with nodes colored by secondary structure.
        
        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to visualize
        pdb_id : str, optional
            PDB identifier for title
        output_path : str, optional
            Path to save the visualization
        figsize : tuple, optional
            Figure size as (width, height)
            
        Returns:
        --------
        None
            Displays or saves the plot
        """
        if graph.number_of_nodes() == 0:
            print("Cannot visualize empty graph!")
            return

        # Set up figure
        plt.figure(figsize=figsize)

        # Collect node positions and colors
        pos = {}
        colors = []
        sizes = []

        # Color map for secondary structure
        ss_colors = {
            'H': 'red',       # Alpha helix
            'B': 'blue',      # Beta bridge
            'E': 'cyan',      # Extended strand
            'G': 'orange',    # 3-10 helix
            'I': 'purple',    # Pi helix
            'T': 'yellow',    # Turn
            'S': 'green',     # Bend
            ' ': 'lightgray', # Coil
            '-': 'gray',      # None
            '?': 'gray'       # Unknown
        }

        # Determine if 'coords' attribute exists
        has_coords = False
        for node in list(graph.nodes())[:5]:  # Check first few nodes
            if 'coords' in graph.nodes[node]:
                has_coords = True
                break

        # For each node, get position and color
        for node in graph.nodes():
            # Get position
            if has_coords and 'coords' in graph.nodes[node]:
                coords = graph.nodes[node]['coords']
                if isinstance(coords, np.ndarray) and coords.size >= 2:
                    pos[node] = (coords[0], coords[1])
                elif isinstance(coords, list) and len(coords) >= 2:
                    pos[node] = (coords[0], coords[1])
            elif 'CA' in graph.nodes[node]:
                # Use CA atom position if coords not available
                ca_pos = graph.nodes[node]['CA']
                if isinstance(ca_pos, np.ndarray) and ca_pos.size >= 2:
                    pos[node] = (ca_pos[0], ca_pos[1])
                elif isinstance(ca_pos, list) and len(ca_pos) >= 2:
                    pos[node] = (ca_pos[0], ca_pos[1])

            # Add color based on secondary structure
            ss = graph.nodes[node].get('ss', '?')
            colors.append(ss_colors.get(ss, 'gray'))

            # Node size
            sizes.append(30)  # Default size

        # If we couldn't get positions for nodes, use spring layout
        if len(pos) < graph.number_of_nodes() * 0.5:  # If less than half the nodes have positions
            print("Using spring layout as fallback")
            pos = nx.spring_layout(graph, seed=42)

        # Draw the graph
        nx.draw(
            graph,
            pos=pos,
            node_color=colors,
            node_size=sizes,
            with_labels=False,
            edge_color='gray',
            width=0.5,
            alpha=0.7
        )

        # Add title
        title = f"Protein Graph (Secondary Structure): {pdb_id}" if pdb_id else "Protein Graph (Secondary Structure)"
        plt.title(title)
        plt.axis('off')

        # Add legend for secondary structure
        legend_handles = []
        ss_names = {
            'H': 'Alpha Helix',
            'B': 'Beta Bridge',
            'E': 'Extended Strand',
            'G': '3-10 Helix',
            'I': 'Pi Helix',
            'T': 'Turn',
            'S': 'Bend',
            ' ': 'Coil',
            '-': 'None',
            '?': 'Unknown'
        }

        for ss, color in ss_colors.items():
            if color in colors:  # Only add to legend if the color is used
                legend_handles.append(mlines.Line2D([], [], color=color, marker='o',
                                                    linestyle='None', markersize=10,
                                                    label=f"{ss} ({ss_names.get(ss, ss)})"))

        plt.legend(handles=legend_handles)

        # Save or show
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def inspect_protein_graph(self, graph, pdb_id=None):
        """
        Inspect a protein graph to verify its structure and properties.
        
        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph object to inspect
        pdb_id : str, optional
            PDB ID for display purposes
            
        Returns:
        --------
        dict
            Dictionary with graph inspection results
        """
        inspection_results = {}

        if pdb_id:
            print(f"\n==== INSPECTING PROTEIN GRAPH: {pdb_id} ====\n")
            inspection_results['pdb_id'] = pdb_id
        else:
            print("\n==== INSPECTING PROTEIN GRAPH ====\n")

        # Basic graph information
        print(f"Graph type: {type(graph)}")
        print(f"Number of nodes: {graph.number_of_nodes()}")
        print(f"Number of edges: {graph.number_of_edges()}")

        inspection_results['graph_type'] = str(type(graph))
        inspection_results['num_nodes'] = graph.number_of_nodes()
        inspection_results['num_edges'] = graph.number_of_edges()

        # Check node information
        inspection_results['nodes'] = {}

        if graph.number_of_nodes() > 0:
            print("\n-- NODE INFORMATION --")
            # Get a sample node
            sample_node = next(iter(graph.nodes()))
            print(f"\nSample node ID: {sample_node}")
            print("Node attributes:")

            inspection_results['nodes']['sample_node'] = sample_node
            inspection_results['nodes']['attributes'] = {}

            for key, value in graph.nodes[sample_node].items():
                if isinstance(value, np.ndarray) and value.size > 10:
                    print(f"  {key}: <ndarray shape={value.shape}, dtype={value.dtype}>")
                    inspection_results['nodes']['attributes'][key] = f"<ndarray shape={value.shape}, dtype={value.dtype}>"
                else:
                    print(f"  {key}: {value}")
                    inspection_results['nodes']['attributes'][key] = str(value)

            # Count of nodes by chain
            chains = {}
            for node in graph.nodes():
                chain = graph.nodes[node].get('chain_id', 'Unknown')
                chains[chain] = chains.get(chain, 0) + 1

            print("\nNodes by chain:")
            inspection_results['nodes']['by_chain'] = {}

            for chain, count in chains.items():
                print(f"  Chain {chain}: {count} nodes")
                inspection_results['nodes']['by_chain'][chain] = count

        # Check edge information
        inspection_results['edges'] = {}

        if graph.number_of_edges() > 0:
            print("\n-- EDGE INFORMATION --")
            # Get a sample edge
            sample_edge = next(iter(graph.edges()))
            print(f"\nSample edge: {sample_edge}")
            print("Edge attributes:")

            inspection_results['edges']['sample_edge'] = str(sample_edge)
            inspection_results['edges']['attributes'] = {}

            for key, value in graph.edges[sample_edge].items():
                print(f"  {key}: {value}")
                inspection_results['edges']['attributes'][key] = str(value)

            # Count edge types
            edge_types = {}
            for u, v, data in graph.edges(data=True):
                # Handle different possible formats of 'kind' attribute
                kind = data.get('kind', 'Unknown')

                if isinstance(kind, set):
                    # If kind is a set, convert to tuple for hashing
                    for k in kind:
                        edge_types[k] = edge_types.get(k, 0) + 1
                elif isinstance(kind, list):
                    # If kind is a list, count each type
                    for k in kind:
                        edge_types[k] = edge_types.get(k, 0) + 1
                else:
                    # Otherwise use as-is
                    edge_types[kind] = edge_types.get(kind, 0) + 1

            print("\nEdge types:")
            inspection_results['edges']['by_type'] = {}

            for edge_type, count in edge_types.items():
                print(f"  {edge_type}: {count} edges")
                inspection_results['edges']['by_type'][str(edge_type)] = count

        # Check for secondary structure
        ss_present = 0
        ss_types = {}
        for node, data in graph.nodes(data=True):
            if 'ss' in data and data['ss'] != '?':
                ss_present += 1
                ss_types[data['ss']] = ss_types.get(data['ss'], 0) + 1

        inspection_results['secondary_structure'] = {
            'present': ss_present > 0,
            'count': ss_present,
            'types': {}
        }

        if ss_present > 0:
            print(f"\nSecondary structure information: Present in {ss_present} nodes")
            print("Secondary structure types:")

            ss_names = {
                'H': 'Alpha Helix',
                'B': 'Beta Bridge',
                'E': 'Extended Strand',
                'G': '3-10 Helix',
                'I': 'Pi Helix',
                'T': 'Turn',
                'S': 'Bend',
                ' ': 'Coil',
                '-': 'None'
            }

            for ss_type, count in ss_types.items():
                name = ss_names.get(ss_type, ss_type)
                print(f"  {ss_type} ({name}): {count} nodes")
                inspection_results['secondary_structure']['types'][ss_type] = {
                    'name': name,
                    'count': count
                }
        else:
            print(f"\nSecondary structure information: Missing")

        # Check for backbone atoms
        backbone_present = 0
        for node, data in graph.nodes(data=True):
            if data.get('has_backbone', False):
                backbone_present += 1

        inspection_results['backbone'] = {
            'present': backbone_present > 0,
            'count': backbone_present
        }

        if backbone_present > 0:
            print(f"Backbone atom information: Present in {backbone_present} nodes")
        else:
            print(f"Backbone atom information: Missing")

        print("\n==== INSPECTION COMPLETE ====\n")

        return inspection_results

    def generate_inspection_report(self, graph, pdb_id=None, output_path=None):
        """
        Generate a detailed inspection report for a protein graph.
        
        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to inspect
        pdb_id : str, optional
            PDB identifier
        output_path : str, optional
            Path to save the report
            
        Returns:
        --------
        dict
            Dictionary with inspection results
        """
        # Perform inspection
        results = self.inspect_protein_graph(graph, pdb_id)

        # Save report if output path is provided
        if output_path:
            import json
            try:
                with open(output_path, 'w') as f:
                    json.dump(results, f, indent=2)
                print(f"Saved inspection report to {output_path}")
            except Exception as e:
                print(f"Error saving inspection report: {str(e)}")

        return results

    def test_structure_info(self, pdb_path, protein_preprocessor):
        """
        Test secondary structure and backbone information extraction.
        
        Parameters:
        -----------
        pdb_path : str
            Path to PDB file
        protein_preprocessor : ProteinPreprocessor
            Preprocessor instance to use for parsing
            
        Returns:
        --------
        networkx.Graph
            Graph with structure information
        """
        print(f"Testing secondary structure and backbone info for: {pdb_path}")
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]

        # Parse PDB
        atoms, residues_by_chain = protein_preprocessor.parse_pdb(pdb_path)
        print(f"Found {len(atoms)} atoms across {len(residues_by_chain)} chains")

        # Calculate secondary structure
        ss_data = protein_preprocessor.calculate_secondary_structure(pdb_path)
        print(f"Found secondary structure for {len(ss_data)} residues")

        # Create graph
        nx_graph = protein_preprocessor.create_protein_graph(pdb_path)
        print(f"Created graph with {nx_graph.number_of_nodes()} nodes and {nx_graph.number_of_edges()} edges")

        # Add structure info
        nx_graph = protein_preprocessor.add_structure_to_graph(nx_graph, ss_data)

        # If backbone info is still missing, try direct extraction from PDB
        backbone_present = 0
        for node, data in nx_graph.nodes(data=True):
            if data.get('has_backbone', False):
                backbone_present += 1

        print(f"Nodes with backbone info after first method: {backbone_present}")

        if backbone_present == 0:
            print("Using PDB file to extract backbone information...")
            nx_graph = protein_preprocessor.add_backbone_from_pdb(nx_graph, pdb_path)

            backbone_present = 0
            for node, data in nx_graph.nodes(data=True):
                if data.get('has_backbone', False):
                    backbone_present += 1

            print(f"Nodes with backbone info after PDB extraction: {backbone_present}")

        # Count SS info
        ss_present = 0
        for node, data in nx_graph.nodes(data=True):
            if 'ss' in data and data['ss'] != '?':
                ss_present += 1

        print(f"Nodes with SS info: {ss_present}")

        # Inspect graph
        self.inspect_protein_graph(nx_graph, pdb_id)

        return nx_graph