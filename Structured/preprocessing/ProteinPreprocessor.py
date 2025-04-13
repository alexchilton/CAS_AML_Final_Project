import os
import numpy as np
import networkx as nx
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import pickle
import json
from tqdm import tqdm
from graphein.protein.config import ProteinGraphConfig, DSSPConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.features.nodes.dssp import add_dssp_df
from graphein.protein.utils import download_pdb
import graphein.protein as gp
from functools import partial
from graphein.ml.conversion import GraphFormatConvertor
from graphein.protein.edges.distance import (add_peptide_bonds,
                                             add_hydrogen_bond_interactions,
                                             add_disulfide_interactions,
                                             add_ionic_interactions,
                                             add_aromatic_interactions,
                                             add_aromatic_sulphur_interactions,
                                             add_cation_pi_interactions
                                             )


class ProteinPreprocessor:
    """
    A class for preprocessing protein structures from PDB files.

    This class handles:
    - PDB parsing
    - Computing relative positions
    - Extracting bond information (bond length, angles, torsions)
    - Assigning partial charges
    - Generating secondary structure labels
    - Graph construction
    """

    def __init__(self, output_dir=None):
        """
        Initialize the protein preprocessor.

        Parameters:
        -----------
        output_dir : str, optional
            Directory to save processed data. If None, will use a 'processed_data'
            subdirectory in the input folder.
        """
        self.output_dir = output_dir
        # Use this edge function set for complete biochemical interactions
        self.all_edge_func = {"edge_construction_functions": [

            add_peptide_bonds,
            add_aromatic_interactions,
            add_hydrogen_bond_interactions,
            add_disulfide_interactions,
            add_ionic_interactions,
            add_aromatic_sulphur_interactions,
            add_cation_pi_interactions,
            gp.add_hydrophobic_interactions,
            gp.add_salt_bridges]}
        # not sure the salt bridges always work


        # Use these metadata configurations for the properties you need
        self.complete_config = {
            #"graph_metadata_functions": [gp.rsa, gp.secondary_structure], - these come from the dssp
            "node_metadata_functions": [gp.amino_acid_one_hot,
                                gp.meiler_embedding,
                                partial(gp.expasy_protein_scale, add_separate=True)]
        }

        # Combined configuration
        self.full_config = gp.ProteinGraphConfig(**{**self.all_edge_func, **self.complete_config})


    def process_pdb_folder(self, folder_path, output_path=None):
        """
        Process all PDB files in a folder and prepare data for PINN GNN.

        Parameters:
        -----------
        folder_path : str
            Path to the folder containing PDB files
        output_path : str, optional
            Optional path to save processed data. If None, will use
            folder_path/processed_data

        Returns:
        --------
        dict
            Dictionary containing processed data for all PDB files
        """
        if output_path is None:
            output_path = os.path.join(folder_path, "processed_data")

        # Create output directory if it doesn't exist
        os.makedirs(output_path, exist_ok=True)

        # Get all PDB files in the folder
        pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]

        # Dictionary to store all processed data
        all_data = {}

        # Process each PDB file
        for pdb_file in tqdm(pdb_files, desc="Processing PDB files"):
            pdb_path = os.path.join(folder_path, pdb_file)
            pdb_id = os.path.splitext(pdb_file)[0]  # PDB ID without extension

            try:
                # Parse PDB file
                atoms, residues_by_chain = self.parse_pdb(pdb_path)

                # Calculate secondary structure
                ss_data = self.calculate_secondary_structure(pdb_path)

                # Create protein graph
                nx_graph = self.create_protein_graph(pdb_path)

                # Add secondary structure and backbone information
                nx_graph = self.add_structure_to_graph(nx_graph, ss_data)

                # Create PDB data with graph
                pdb_data = {
                    'pdb_id': pdb_id,
                    'atoms': atoms,
                    'residues_by_chain': residues_by_chain,
                    'nx_graph': nx_graph,
                    'relative_positions': {},
                    'edge_indices': {},
                    'bond_lengths': {},
                    'angles': {},
                    'torsions': {},
                    'secondary_structure': {},
                    'backbone_atoms': {},
                    'charges': {},
                    'hydrophobic_residues': {}
                }

                # Process each chain
                for chain_id, chain_residues in residues_by_chain.items():
                    # Get atoms for this chain
                    chain_atoms = [atom for atom in atoms if atom[0] == chain_id]

                    # Skip empty chains
                    if not chain_atoms:
                        continue

                    # Extract positions
                    positions = np.array([atom[5] for atom in chain_atoms])

                    # Compute relative positions
                    rel_pos = self.compute_relative_positions_chain(positions)

                    # Create edge index for graph construction
                    edge_index = self.create_edge_index(chain_atoms)

                    # Get backbone atoms
                    backbone_atoms = self.extract_backbone_atoms(chain_atoms)

                    # Get secondary structure for this chain
                    chain_ss = {res_id: ss_data.get((chain_id, res_id), '?')
                                for res_id, _ in chain_residues}

                    # Calculate physical properties
                    bond_lengths = self.calculate_bond_lengths(chain_atoms, edge_index)
                    angles = self.calculate_angles(chain_atoms, edge_index)
                    torsions = self.calculate_torsions(chain_atoms, edge_index)
                    charges = self.assign_charges(chain_atoms)
                    hydrophobic = self.identify_hydrophobic_residues(chain_residues)

                    # Store in PDB data dictionary
                    pdb_data['relative_positions'][chain_id] = rel_pos
                    pdb_data['edge_indices'][chain_id] = edge_index
                    pdb_data['bond_lengths'][chain_id] = bond_lengths
                    pdb_data['angles'][chain_id] = angles
                    pdb_data['torsions'][chain_id] = torsions
                    pdb_data['secondary_structure'][chain_id] = chain_ss
                    pdb_data['backbone_atoms'][chain_id] = backbone_atoms
                    pdb_data['charges'][chain_id] = charges
                    pdb_data['hydrophobic_residues'][chain_id] = hydrophobic

                # Store in the main dictionary
                all_data[pdb_id] = pdb_data

                # Save processed data
                if output_path:
                    self.save_processed_data(pdb_data, pdb_id, output_path)

            except Exception as e:
                print(f"Error processing {pdb_file}: {str(e)}")
                import traceback
                traceback.print_exc()

        # Save summary of all processed data
        self.save_processing_summary(all_data, output_path)

        return all_data

    def save_processed_data(self, pdb_data, pdb_id, output_path):
        """
        Save processed data to disk.

        Parameters:
        -----------
        pdb_data : dict
            Processed data for a single PDB
        pdb_id : str
            PDB identifier
        output_path : str
            Directory to save the data
        """
        # Save the NetworkX graph
        graph_path = os.path.join(output_path, f"{pdb_id}_graph.pkl")
        try:
            with open(graph_path, 'wb') as f:
                pickle.dump(pdb_data['nx_graph'], f)
            print(f"Saved graph to {graph_path}")
        except Exception as e:
            print(f"Error saving graph for {pdb_id}: {str(e)}")

        # Save other data as numpy file
        # Remove graph from data to save separately (graphs can be large)
        data_to_save = {k: v for k, v in pdb_data.items() if k != 'nx_graph'}
        data_path = os.path.join(output_path, f"{pdb_id}_data.npy")
        try:
            np.save(data_path, data_to_save)
            print(f"Saved data to {data_path}")
        except Exception as e:
            print(f"Error saving data for {pdb_id}: {str(e)}")

    def save_processing_summary(self, all_data, output_path):
        """
        Save a summary of the processing results.

        Parameters:
        -----------
        all_data : dict
            Dictionary containing all processed PDB data
        output_path : str
            Directory to save the summary
        """
        summary = {}
        for pdb_id, pdb_data in all_data.items():
            summary[pdb_id] = {
                'num_atoms': len(pdb_data['atoms']),
                'num_chains': len(pdb_data['residues_by_chain']),
                'chains': list(pdb_data['residues_by_chain'].keys()),
                'num_nodes': pdb_data['nx_graph'].number_of_nodes(),
                'num_edges': pdb_data['nx_graph'].number_of_edges(),
                'has_secondary_structure': any(len(pdb_data['secondary_structure'].get(chain, {})) > 0
                                               for chain in pdb_data['residues_by_chain']),
                'has_backbone_atoms': any(len(pdb_data['backbone_atoms'].get(chain, {})) > 0
                                          for chain in pdb_data['residues_by_chain'])
            }

        summary_path = os.path.join(output_path, "processing_summary.json")
        try:
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2)
            print(f"Saved processing summary to {summary_path}")
        except Exception as e:
            print(f"Error saving summary: {str(e)}")

    def process_single_pdb(self, pdb_path):
        """
        Process a single PDB file.

        Parameters:
        -----------
        pdb_path : str
            Path to the PDB file

        Returns:
        --------
        dict
            Dictionary containing processed data
        """
        pdb_id = os.path.splitext(os.path.basename(pdb_path))[0]

        # Parse PDB file
        atoms, residues_by_chain = self.parse_pdb(pdb_path)

        # Calculate secondary structure
        ss_data = self.calculate_secondary_structure(pdb_path)

        # Create protein graph
        nx_graph = self.create_protein_graph(pdb_path)

        # Add secondary structure and backbone information
        nx_graph = self.add_structure_to_graph(nx_graph, ss_data)

        # Create PDB data with graph
        pdb_data = {
            'pdb_id': pdb_id,
            'atoms': atoms,
            'residues_by_chain': residues_by_chain,
            'nx_graph': nx_graph,
            'relative_positions': {},
            'edge_indices': {},
            'bond_lengths': {},
            'angles': {},
            'torsions': {},
            'secondary_structure': {},
            'backbone_atoms': {},
            'charges': {},
            'hydrophobic_residues': {}
        }

        # Process each chain
        for chain_id, chain_residues in residues_by_chain.items():
            # Get atoms for this chain
            chain_atoms = [atom for atom in atoms if atom[0] == chain_id]

            # Skip empty chains
            if not chain_atoms:
                continue

            # Extract positions
            positions = np.array([atom[5] for atom in chain_atoms])

            # Compute relative positions
            rel_pos = self.compute_relative_positions_chain(positions)

            # Create edge index for graph construction
            edge_index = self.create_edge_index(chain_atoms)

            # Get backbone atoms
            backbone_atoms = self.extract_backbone_atoms(chain_atoms)

            # Get secondary structure for this chain
            chain_ss = {res_id: ss_data.get((chain_id, res_id), '?')
                        for res_id, _ in chain_residues}

            # Calculate physical properties
            bond_lengths = self.calculate_bond_lengths(chain_atoms, edge_index)
            angles = self.calculate_angles(chain_atoms, edge_index)
            torsions = self.calculate_torsions(chain_atoms, edge_index)
            charges = self.assign_charges(chain_atoms)
            hydrophobic = self.identify_hydrophobic_residues(chain_residues)

            # Store in PDB data dictionary
            pdb_data['relative_positions'][chain_id] = rel_pos
            pdb_data['edge_indices'][chain_id] = edge_index
            pdb_data['bond_lengths'][chain_id] = bond_lengths
            pdb_data['angles'][chain_id] = angles
            pdb_data['torsions'][chain_id] = torsions
            pdb_data['secondary_structure'][chain_id] = chain_ss
            pdb_data['backbone_atoms'][chain_id] = backbone_atoms
            pdb_data['charges'][chain_id] = charges
            pdb_data['hydrophobic_residues'][chain_id] = hydrophobic

        return pdb_data

    def parse_pdb(self, pdb_filename):
        """
        Parse a PDB file to extract atom and residue information.

        Parameters:
        -----------
        pdb_filename : str
            Path to the PDB file

        Returns:
        --------
        tuple
            (atoms, residues_by_chain) where:
            - atoms is a list of (chain_id, residue_id, residue_name, atom_name, element, position)
            - residues_by_chain is a dict mapping chain IDs to lists of (residue_id, residue_name)
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("Protein", pdb_filename)

        atoms = []
        residues = []

        # Keep track of processed residues by chain
        all_residues_by_chain = {}

        for model in structure:
            for chain in model:
                chain_id = chain.get_id()
                chain_residues = []

                for residue in chain:
                    # Get residue details
                    residue_id = residue.get_id()[1]  # The sequence number
                    residue_name = residue.get_resname()

                    # Add residue to this chain's list
                    chain_residues.append((residue_id, residue_name))

                    # Process all atoms in this residue
                    for atom in residue:
                        atom_name = atom.get_name()
                        element = atom.element
                        position = atom.get_coord()
                        # Include chain ID in atom information for clarity
                        atoms.append((chain_id, residue_id, residue_name, atom_name, element, position))

                # Store all residues for this chain
                all_residues_by_chain[chain_id] = chain_residues

        return atoms, all_residues_by_chain

    def calculate_secondary_structure(self, pdb_path):
        """
        Calculate secondary structure using DSSP.

        Parameters:
        -----------
        pdb_path : str
            Path to the PDB file

        Returns:
        --------
        dict
            Dictionary mapping (chain_id, residue_id) to secondary structure
        """
        try:
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure("protein", pdb_path)
            model = structure[0]

            try:
                # Try direct DSSP call
                dssp = DSSP(model, pdb_path)
            except Exception as e1:
                try:
                    # Fallback: Try with DSSP executable if installed
                    from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                    dssp_dict, dssp_keys = dssp_dict_from_pdb_file(pdb_path)

                    # Convert to format similar to DSSP object
                    dssp = {}
                    for key in dssp_keys:
                        dssp[key] = dssp_dict[key]
                except Exception as e2:
                    print(f"DSSP calculation failed: {str(e1)} | {str(e2)}")
                    # If all methods fail, use a simple fallback
                    return {}

            ss_data = {}

            for key in dssp.keys():
                if isinstance(key, tuple):
                    # Newer DSSP version
                    chain_id = key[0]
                    res_id = key[1][1]  # Residue number
                    if isinstance(dssp[key], dict):
                        ss = dssp[key].get('secstruct', '?')
                    else:
                        # Older DSSP version
                        ss = dssp[key][2]  # Secondary structure
                else:
                    # Alternative key format
                    chain_id, res_id = key[0], key[1]
                    ss = dssp[key][2]

                ss_data[(chain_id, res_id)] = ss

            return ss_data
        except Exception as e:
            print(f"Secondary structure calculation failed: {str(e)}")
            return {}

    def create_protein_graph(self, pdb_path):
        """
        Create a protein graph using Graphein or fallback to custom implementation.

        Parameters:
        -----------
        pdb_path : str
            Path to the PDB file

        Returns:
        --------
        networkx.Graph
            NetworkX graph representation of the protein
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

            # Add edges between consecutive residues and close residues
            # First, sort nodes by chain and residue number
            sorted_nodes = sorted(graph.nodes(),
                                  key=lambda x: (x.split(':')[0],
                                                 int(x.split(':')[2]) if x.split(':')[2].isdigit() else 0))

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

    def add_structure_to_graph(self, graph, ss_data):
        """
        Add secondary structure and backbone information to a protein graph.

        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to update
        ss_data : dict
            Dictionary of secondary structure data

        Returns:
        --------
        networkx.Graph
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

            # Add backbone information
            graph.nodes[node]['has_backbone'] = backbone_atoms_present

        return graph

    def add_backbone_from_pdb(self, graph, pdb_path):
        """
        Add backbone information to the graph by directly extracting it from the PDB file.

        Parameters:
        -----------
        graph : networkx.Graph
            NetworkX graph to update
        pdb_path : str
            Path to the PDB file

        Returns:
        --------
        networkx.Graph
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

    def extract_backbone_atoms(self, chain_atoms):
        """
        Extract backbone atoms (N, CA, C, O) from a chain.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain

        Returns:
        --------
        dict
            Dictionary mapping residue IDs to backbone atom positions
        """
        backbone = {}

        for atom in chain_atoms:
            chain_id, res_id, res_name, atom_name, element, position = atom

            if atom_name in ['N', 'CA', 'C', 'O']:
                if res_id not in backbone:
                    backbone[res_id] = {}
                backbone[res_id][atom_name] = position

        return backbone

    def compute_relative_positions_chain(self, positions):
        """
        Compute relative positions between atoms in a chain.

        Parameters:
        -----------
        positions : numpy.ndarray
            Array of atom positions

        Returns:
        --------
        numpy.ndarray
            Array of relative positions
        """
        n_atoms = positions.shape[0]
        rel_pos = np.zeros((n_atoms, n_atoms, 3))

        for i in range(n_atoms):
            for j in range(n_atoms):
                rel_pos[i, j] = positions[j] - positions[i]

        return rel_pos

    def create_edge_index(self, chain_atoms):
        """
        Create edge index for graph construction.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain

        Returns:
        --------
        numpy.ndarray
            Array of edge indices
        """
        n_atoms = len(chain_atoms)
        edges = []

        # Example: Connect atoms within the same residue and consecutive residues
        for i in range(n_atoms):
            res_i = chain_atoms[i][1]  # Residue ID
            for j in range(i+1, n_atoms):
                res_j = chain_atoms[j][1]  # Residue ID
                if res_i == res_j or abs(res_i - res_j) == 1:
                    edges.append([i, j])
                    edges.append([j, i])  # For undirected graph

        return np.array(edges).T if edges else np.zeros((2, 0), dtype=int)

    def calculate_bond_lengths(self, chain_atoms, edge_index):
        """
        Calculate bond lengths between connected atoms.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain
        edge_index : numpy.ndarray
            Array of edge indices

        Returns:
        --------
        dict
            Dictionary of bond lengths
        """
        # Placeholder for bond length calculation
        # You would implement actual bond length calculation here
        return {}

    def calculate_angles(self, chain_atoms, edge_index):
        """
        Calculate angles between connected atoms.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain
        edge_index : numpy.ndarray
            Array of edge indices

        Returns:
        --------
        dict
            Dictionary of angles
        """
        # Placeholder for angle calculation
        # You would implement actual angle calculation here
        return {}

    def calculate_torsions(self, chain_atoms, edge_index):
        """
        Calculate torsion angles between connected atoms.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain
        edge_index : numpy.ndarray
            Array of edge indices

        Returns:
        --------
        dict
            Dictionary of torsion angles
        """
        # Placeholder for torsion calculation
        # You would implement actual torsion calculation here
        return {}

    def assign_charges(self, chain_atoms):
        """
        Assign partial charges to atoms based on residue type and atom name.

        Parameters:
        -----------
        chain_atoms : list
            List of atoms in a chain

        Returns:
        --------
        dict
            Dictionary of partial charges
        """
        charges = {}

        # Simplified charge assignment based on atom types
        charge_map = {
            'N': -0.5,   # Backbone nitrogen
            'C': 0.5,    # Backbone carbon
            'O': -0.5,   # Carbonyl oxygen
            'S': 0.0,    # Sulfur
            'P': 0.5     # Phosphorus
        }

        for i, atom in enumerate(chain_atoms):
            chain_id, res_id, res_name, atom_name, element, position = atom
            charges[i] = charge_map.get(element, 0.0)

        return charges

    def identify_hydrophobic_residues(self, chain_residues):
        """
        Identify hydrophobic residues in a chain.

        Parameters:
        -----------
        chain_residues : list
            List of residues in a chain

        Returns:
        --------
        set
            Set of residue IDs that are hydrophobic
        """
        hydrophobic_aas = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR'}
        hydrophobic_res_ids = set()

        for res_id, res_name in chain_residues:
            if res_name in hydrophobic_aas:
                hydrophobic_res_ids.add(res_id)

        return hydrophobic_res_ids


    def generate_graph_no_dssp(self, path = "2X89.pdb"):
        """
        Generate a graph without using DSSP.

        Parameters:
        -----------
        path : str
            Path to the PDB file

        Returns:
        --------
        networkx.Graph
            NetworkX graph representation of the protein
        """

        graph = construct_graph(path=path, config=self.full_config)
        return graph