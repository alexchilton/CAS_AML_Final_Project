import numpy as np
from Bio import PDB
import networkx as nx

def create_protein_graph(pdb_path, ss_data=None):
    """
    Create a protein graph using Graphein or fallback to custom implementation.
    
    Parameters:
    pdb_path -- Path to the PDB file
    ss_data -- Optional secondary structure data dictionary
    
    Returns:
    NetworkX graph
    """
    try:
        # Try the current Graphein API
        import graphein.protein as gp
        
        # Check the actual function signature
        import inspect
        sig = inspect.signature(gp.construct_graph)
        
        # Check if 'pdb_path' is a parameter
        if 'pdb_path' in sig.parameters:
            graph = gp.construct_graph(pdb_path=pdb_path)
        elif 'path' in sig.parameters:
            # Try with 'path' instead
            graph = gp.construct_graph(path=pdb_path)
        else:
            # Check if the function expects the first argument to be the path
            params = list(sig.parameters.keys())
            if len(params) > 0:
                # Try with positional argument
                graph = gp.construct_graph(pdb_path)
            else:
                raise ValueError("Cannot determine correct parameters for construct_graph")
                
        return graph
        
    except Exception as e:
        print(f"Graphein error: {str(e)}")
        print("Falling back to custom graph implementation...")
        
        # Create a simple graph using NetworkX and BioPython
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        
        # Create basic graph
        graph = nx.Graph()
        
        # Add nodes for each residue
        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    # Skip non-standard residues
                    if residue.get_id()[0] != " ":
                        continue
                        
                    res_id = residue.get_id()[1]
                    res_name = residue.get_resname()
                    
                    # Add node
                    node_id = f"{chain_id}:{res_name}:{res_id}"
                    graph.add_node(
                        node_id,
                        chain_id=chain_id,
                        residue_number=res_id,
                        residue_name=res_name
                    )
                    
                    # Add atom data and track if backbone atoms are present
                    has_backbone = False
                    for atom in residue:
                        atom_name = atom.get_name()
                        coord = atom.get_coord()
                        
                        # Store coordinates of backbone atoms
                        if atom_name in ['CA', 'N', 'C', 'O']:
                            graph.nodes[node_id][atom_name] = coord
                            has_backbone = True
                    
                    # Set backbone flag
                    graph.nodes[node_id]['has_backbone'] = has_backbone
                    
                    # Store CA coordinates for distance calculations
                    if 'CA' in graph.nodes[node_id]:
                        graph.nodes[node_id]['coords'] = graph.nodes[node_id]['CA']
                    
                    # Add secondary structure if available
                    if ss_data and (chain_id, res_id) in ss_data:
                        graph.nodes[node_id]['ss'] = ss_data[(chain_id, res_id)]
                    else:
                        graph.nodes[node_id]['ss'] = '?'
        
        # Add edges between consecutive residues and close residues
        # First, sort nodes by chain and residue number
        sorted_nodes = sorted(graph.nodes(), key=lambda x: (x.split(':')[0], int(x.split(':')[2]) if x.split(':')[2].isdigit() else 0))
        
        # Add backbone connections
        for i in range(len(sorted_nodes) - 1):
            node1 = sorted_nodes[i]
            chain1, _, res1 = node1.split(':')
            
            for j in range(i + 1, min(i + 20, len(sorted_nodes))):  # Check only nearby residues
                node2 = sorted_nodes[j]
                chain2, _, res2 = node2.split(':')
                
                # Add edge if same chain and sequential
                if chain1 == chain2 and abs(int(res1) - int(res2)) == 1:
                    graph.add_edge(node1, node2, kind={'peptide_bond'})
                
                # Add edge if CA atoms are within 7 Angstroms
                if 'coords' in graph.nodes[node1] and 'coords' in graph.nodes[node2]:
                    ca1 = graph.nodes[node1]['coords']
                    ca2 = graph.nodes[node2]['coords']
                    distance = np.linalg.norm(ca1 - ca2)
                    
                    if distance < 7.0:  # 7 Angstrom cutoff
                        graph.add_edge(node1, node2, kind={'contact'}, distance=distance)
        
        return graph

def add_structure_to_graph(graph, ss_data):
    """
    Add secondary structure and backbone information to a protein graph.
    
    Parameters:
    graph -- NetworkX graph to update
    ss_data -- Dictionary of secondary structure data
    
    Returns:
    Updated NetworkX graph
    """
    for node, data in graph.nodes(data=True):
        # Extract chain and residue information
        if ':' in node:  # Format like 'H:VAL:1'
            parts = node.split(':')
            chain_id = parts[0]
            if len(parts) >= 3:
                try:
                    residue_id = int(parts[2])
                except ValueError:
                    residue_id = parts[2]
        else:
            # Try to extract from attributes
            chain_id = data.get('chain_id', '')
            residue_id = data.get('residue_number', None)
        
        # Add secondary structure if available
        if (chain_id, residue_id) in ss_data:
            graph.nodes[node]['ss'] = ss_data[(chain_id, residue_id)]
        else:
            graph.nodes[node]['ss'] = '?'
        
        # Check for backbone atoms using exact attribute names from Graphein
        backbone_keys = ['CA', 'C', 'N', 'O', 'atom_type', 'element_symbol', 'coords']
        backbone_atoms_present = False
        
        for key in data.keys():
            # Check for backbone atom coordinates or any CA atom indicator
            if key in ['CA', 'C', 'N', 'O']:
                backbone_atoms_present = True
                break
            elif key == 'atom_type' and data[key] in ['CA', 'C', 'N', 'O']:
                backbone_atoms_present = True
                break
        
        # For debugging, print the first node's data
        if node == list(graph.nodes())[0]:
            print(f"First node data keys: {list(data.keys())}")
            print(f"Checking for backbone atoms: {backbone_atoms_present}")
            
        # Add backbone information
        graph.nodes[node]['has_backbone'] = backbone_atoms_present
    
    return graph

def add_backbone_from_pdb(graph, pdb_path):
    """
    Add backbone information to the graph by directly extracting it from the PDB file.
    
    Parameters:
    graph -- NetworkX graph to update
    pdb_path -- Path to the PDB file
    
    Returns:
    Updated NetworkX graph
    """
    # Parse PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    # Extract backbone atoms
    backbone_info = {}
    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                if residue.get_id()[0] != " ":  # Skip non-standard residues
                    continue
                    
                res_id = residue.get_id()[1]
                
                # Check for backbone atoms
                backbone_atoms = {}
                for atom_name in ['CA', 'N', 'C', 'O']:
                    if atom_name in residue:
                        backbone_atoms[atom_name] = residue[atom_name].get_coord()
                
                if backbone_atoms:
                    backbone_info[(chain_id, res_id)] = backbone_atoms
    
    # Add backbone info to graph
    for node, data in graph.nodes(data=True):
        # Extract chain and residue ID
        if ':' in node:
            parts = node.split(':')
            chain_id = parts[0]
            if len(parts) >= 3:
                try:
                    res_id = int(parts[2])
                except ValueError:
                    res_id = parts[2]
                    
                # Check if we have backbone info for this residue
                if (chain_id, res_id) in backbone_info:
                    # Add backbone atom coordinates
                    for atom_name, coord in backbone_info[(chain_id, res_id)].items():
                        graph.nodes[node][atom_name] = coord
                    
                    # Set backbone flag
                    graph.nodes[node]['has_backbone'] = True
                else:
                    graph.nodes[node]['has_backbone'] = False
        else:
            # Extract from attributes if possible
            chain_id = data.get('chain_id', '')
            res_id = data.get('residue_number', None)
            
            if chain_id and res_id and (chain_id, res_id) in backbone_info:
                # Add backbone atom coordinates
                for atom_name, coord in backbone_info[(chain_id, res_id)].items():
                    graph.nodes[node][atom_name] = coord
                
                # Set backbone flag
                graph.nodes[node]['has_backbone'] = True
    
    return graph
