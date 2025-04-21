from GraphBuilder import GraphBuilder
import matplotlib.pyplot as plt
import networkx as nx

# Example usage of GraphBuilder class
def debug_graph(graph):
    ''' Debugging function to print basic graph info '''
    # Check graph type
    print(f"Graph type: {type(graph)}")

    # Print basic node info
    print(f"\nNumber of nodes: {len(graph.nodes())}")
    if len(graph.nodes()) > 0:
        sample_node = list(graph.nodes())[0]
        print(f"Sample node: {sample_node}")
        print(f"Sample node attributes: {graph.nodes[sample_node]}")

    # Print basic edge info
    print(f"\nNumber of edges: {len(graph.edges())}")
    if len(graph.edges()) > 0:
        sample_edge = list(graph.edges())[0]
        print(f"Sample edge: {sample_edge}")
        print(f"Sample edge attributes: {graph.edges[sample_edge]}")


def basic_usage_example():
    """Basic example of using GraphBuilder to load and explore graph data"""
    # Initialize GraphBuilder with path to data directory
    builder = GraphBuilder("/Users/alexchilton/Downloads/nanos_networkx_small")

    # List available nanobody IDs
    nanobody_ids = builder.get_nanobody_ids()
    print(f"Found {len(nanobody_ids)} nanobodies")
    print(f"First 5 nanobodies: {nanobody_ids[:5]}")

    # Load NetworkX graph for a specific nanobody
    nanobody_id = nanobody_ids[0]
    graph = builder.load_graph(nanobody_id)

    # Print graph summary
    builder.print_graph_summary(graph)

    debug_graph(graph)
    # List available chain files
    chain_files = builder.list_chain_files(nanobody_id)
    print(f"\nAvailable chain files for {nanobody_id}:")
    for chain_id, files in chain_files.items():
        print(f"  Chain {chain_id}: {len(files)} files")

    # Get available file types for a specific chain
    if chain_files:
        chain_id = list(chain_files.keys())[0]
        file_types = builder.get_available_file_types(nanobody_id, chain_id)
        print(f"\nAvailable file types for chain {chain_id}:")
        print(f"  {file_types}")

def load_chain_data_example():
    """Example of loading specific chain data"""
    builder = GraphBuilder("/Users/alexchilton/Downloads/nanos_networkx_small")
    nanobody_ids = builder.get_nanobody_ids()
    nanobody_id = nanobody_ids[0]

    # Get chain IDs
    chain_files = builder.list_chain_files(nanobody_id)
    if not chain_files:
        print("No chain files found")
        return

    chain_id = list(chain_files.keys())[0]

    # Load secondary structure data
    try:
        ss_data = builder.load_chain_file(nanobody_id, chain_id, "ss")
        print(f"Secondary structure data for {nanobody_id}, chain {chain_id}:")
        print(f"  Number of residues: {len(ss_data)}")
        print(f"  SS types: {set(ss_data.values())}")
    except FileNotFoundError:
        print("Secondary structure file not found")

    # Load backbone data
    try:
        backbone_data = builder.load_chain_file(nanobody_id, chain_id, "backbone")
        print(f"\nBackbone data for {nanobody_id}, chain {chain_id}:")
        print(f"  Number of residues: {len(backbone_data)}")
        if backbone_data:
            # Show sample residue
            res_id = list(backbone_data.keys())[0]
            print(f"  Sample residue {res_id} backbone atoms:")
            for atom, coords in backbone_data[res_id].items():
                print(f"    {atom}: {coords}")
    except FileNotFoundError:
        print("Backbone file not found")

def convert_to_pytorch_geometric_example():
    """Example of converting NetworkX graph to PyTorch Geometric format"""
    builder = GraphBuilder("/Users/alexchilton/Downloads/nanos_networkx_small")
    nanobody_ids = builder.get_nanobody_ids()
    nanobody_id = nanobody_ids[0]

    # Load graph
    graph = builder.load_graph(nanobody_id)

    # Convert to PyTorch Geometric
    data = builder.convert_to_pytorch_geometric(graph)

    print(f"\nConverted {nanobody_id} to PyTorch Geometric:")
    print(f"  Node features shape: {data.x.shape}")
    print(f"  Edge index shape: {data.edge_index.shape}")
    print(f"  Edge attributes shape: {data.edge_attr.shape}")

def visualize_graph_example():
    """Example of visualizing a protein graph"""
    builder = GraphBuilder("/Users/alexchilton/Downloads/nanos_networkx_small")
    nanobody_ids = builder.get_nanobody_ids()
    nanobody_id = nanobody_ids[0]

    # Load graph
    graph = builder.load_graph(nanobody_id)

    # Extract 3D coordinates for visualization
    pos = {}
    for node, attrs in graph.nodes(data=True):
        if 'coords' in attrs:
            pos[node] = attrs['coords'][:2]  # Use only x, y for 2D visualization
        elif 'CA' in attrs:
            pos[node] = attrs['CA'][:2]

    # Create a smaller graph if the original is too large
    if len(graph) > 50:
        # Create a subgraph of the first 50 nodes
        nodes = list(graph.nodes())[:50]
        subgraph = graph.subgraph(nodes)
        visualize_graph = subgraph
        pos = {node: pos[node] for node in nodes if node in pos}
    else:
        visualize_graph = graph

    # Set up colors based on chain ID
    colors = []
    for node in visualize_graph.nodes():
        chain_id = visualize_graph.nodes[node].get('chain_id', '')
        # Simple hash function to assign colors
        color_val = hash(chain_id) % 6
        colors.append(color_val)

    # Draw the graph
    plt.figure(figsize=(10, 8))
    nx.draw(visualize_graph, pos=pos, with_labels=False,
            node_color=colors, node_size=50, cmap=plt.cm.tab10)
    plt.title(f"Graph visualization for {nanobody_id}")
    plt.show()

def process_features_example():
    """Example of using ProteinFeatureProcessor with GraphBuilder"""
    from GraphBuilder import GraphBuilder
    from ProteinFeatureProcessor import ProteinFeatureProcessor

    # Initialize GraphBuilder and ProteinFeatureProcessor
    builder = GraphBuilder("/Users/alexchilton/Downloads/nanos_networkx_small")
    processor = ProteinFeatureProcessor()

    # Get available nanobody IDs
    nanobody_ids = builder.get_nanobody_ids()
    print(f"Found {len(nanobody_ids)} nanobodies")

    # Load a graph
    nanobody_id = nanobody_ids[0]
    print(f"Processing {nanobody_id}")
    graph = builder.load_graph(nanobody_id)

    # Print original graph summary
    print("\nOriginal graph summary:")
    builder.print_graph_summary(graph)

    # Process graph to create PyTorch Geometric data object
    data = processor.process_graph_to_pyg(graph)

    # Print PyTorch Geometric data summary
    print("\nProcessed PyTorch Geometric data:")
    print(f"  Node features shape: {data.x.shape}")
    print(f"  Edge index shape: {data.edge_index.shape}")
    print(f"  Edge attributes shape: {data.edge_attr.shape}")

    # Print feature interpretation guide
    print("\nNode feature interpretation:")
    print(f"  Features 0-{len(processor.amino_acids)-1}: Amino acid one-hot encoding")
    print(f"  Features {len(processor.amino_acids)}-{len(processor.amino_acids)+len(processor.ss_types)-1}: Secondary structure one-hot encoding")
    print(f"  Features {len(processor.amino_acids)+len(processor.ss_types)}-{len(processor.amino_acids)+len(processor.ss_types)+2}: 3D coordinates (x,y,z)")
    print(f"  Feature {len(processor.amino_acids)+len(processor.ss_types)+3}: B-factor")
    print(f"  Features {len(processor.amino_acids)+len(processor.ss_types)+4}-{len(processor.amino_acids)+len(processor.ss_types)+10}: Meiler features")

    # Print edge feature interpretation
    print("\nEdge feature interpretation:")
    print("  Feature 0: Distance between residues")
    print("  Feature 1: Edge type (1.0 for peptide bond, 0.0 for spatial contact)")

    return data

if __name__ == "__main__":
    #basic_usage_example()
    #load_chain_data_example()
    #convert_to_pytorch_geometric_example()
    # Uncomment to run visualization example
    # visualize_graph_example()
    process_features_example()