from graphein.protein.config import ProteinGraphConfig, DSSPConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.features.nodes.dssp import add_dssp_df
from graphein.protein.utils import download_pdb
import graphein.protein as gp
from functools import partial
import networkx as nx
from graphein.ml.conversion import GraphFormatConvertor
from graphein.protein.edges.distance import (add_peptide_bonds,
                                             add_hydrogen_bond_interactions,
                                             add_disulfide_interactions,
                                             add_ionic_interactions,
                                             add_aromatic_interactions,
                                             add_aromatic_sulphur_interactions,
                                             add_cation_pi_interactions
                                             )


class FeatureExtractor:
    """Extracts features from protein structures"""
    # Methods for extracting different types of features
    def __init__(self):
        self.config = ProteinGraphConfig()
        self.config.dssp_config = DSSPConfig()

    def add_dssp_df(self, graph, dssp_config= None):
        if dssp_config is None:
            dssp_config = self.config.dssp_config
        add_dssp_df(graph, dssp_config)
        return graph

    def add_all_dssp_features_to_nodes(self, graph):
        """
        Extract all available DSSP features from the dataframe and add them to nodes.
        """
        if "dssp_df" not in graph.graph:
            print("No DSSP dataframe found in graph")
            return graph

        # Get the DSSP dataframe
        dssp_df = graph.graph["dssp_df"]

        # Print available columns in DSSP dataframe
        print(f"Available DSSP features: {list(dssp_df.columns)}")

        # Map of common DSSP column names and their descriptions
        dssp_features = {
            'ss': 'Secondary structure (ss_value)',
            'aa': 'Amino acid (dssp_aa)',
            'acc': 'Absolute solvent accessibility (acc)',
            'phi': 'Phi angle (phi)',
            'psi': 'Psi angle (psi)',
            'dssp_index': 'DSSP residue index (dssp_index)',
            'NH_O_1_relidx': 'First NH-O hydrogen bond relative index (nh_o1_relidx)',
            'NH_O_1_energy': 'First NH-O hydrogen bond energy (nh_o1_energy)',
            'O_NH_1_relidx': 'First O-NH hydrogen bond relative index (o_nh1_relidx)',
            'O_NH_1_energy': 'First O-NH hydrogen bond energy (o_nh1_energy)',
            'NH_O_2_relidx': 'Second NH-O hydrogen bond relative index (nh_o2_relidx)',
            'NH_O_2_energy': 'Second NH-O hydrogen bond energy (nh_o2_energy)',
            'O_NH_2_relidx': 'Second O-NH hydrogen bond relative index (o_nh2_relidx)',
            'O_NH_2_energy': 'Second O-NH hydrogen bond energy (o_nh2_energy)',
        }

        # Variation in column naming across different versions
        alt_column_names = {
            'ss': ['ss', 'SS', 'sec_struc'],
            'aa': ['aa', 'AA', 'amino_acid'],
            'acc': ['acc', 'ACC', 'accessibility'],
            'phi': ['phi', 'PHI'],
            'psi': ['psi', 'PSI'],
            'dssp_index': ['dssp_index', 'id'],
            'NH_O_1_relidx': ['NH_O_1_relidx', 'NH-O_1_relidx', 'NH_O_1_ridx'],
            'NH_O_1_energy': ['NH_O_1_energy', 'NH-O_1_energy'],
            'O_NH_1_relidx': ['O_NH_1_relidx', 'O-NH_1_relidx', 'O_NH_1_ridx'],
            'O_NH_1_energy': ['O_NH_1_energy', 'O-NH_1_energy'],
            'NH_O_2_relidx': ['NH_O_2_relidx', 'NH-O_2_relidx', 'NH_O_2_ridx'],
            'NH_O_2_energy': ['NH_O_2_energy', 'NH-O_2_energy'],
            'O_NH_2_relidx': ['O_NH_2_relidx', 'O-NH_2_relidx', 'O_NH_2_ridx'],
            'O_NH_2_energy': ['O_NH_2_energy', 'O-NH_2_energy'],
        }

        # Find the actual column names in the dataframe
        actual_columns = {}
        for feature, alternatives in alt_column_names.items():
            for alt in alternatives:
                if alt in dssp_df.columns:
                    actual_columns[feature] = alt
                    break

        print(f"Found {len(actual_columns)} DSSP features in the dataframe")

        # Track how many nodes were updated
        updated_nodes = 0
        features_added = set()

        # Different versions of Graphein might have different column names for chain & residue
        # Try to identify the correct column names
        chain_col = next((c for c in dssp_df.columns if 'chain' in c.lower()), None)
        res_num_col = next((c for c in dssp_df.columns if 'res' in c.lower() and 'num' in c.lower()), None)

        if not chain_col or not res_num_col:
            print(f"Could not identify chain and residue number columns in DSSP dataframe")
            print(f"Available columns: {dssp_df.columns}")
            return graph

        print(f"Using '{chain_col}' for chain ID and '{res_num_col}' for residue number")

        # Add features to each node
        for node, data in graph.nodes(data=True):
            # Extract chain and residue info from node ID
            parts = str(node).split(':')
            if len(parts) < 3:
                continue

            chain = parts[0]
            residue_num = parts[2]

            try:
                # Find matching row in DSSP dataframe
                mask = (dssp_df[chain_col] == chain) & (dssp_df[res_num_col] == int(residue_num))
                matching_rows = dssp_df[mask]

                if matching_rows.empty:
                    continue

                # Add all available features to the node
                for feature, col_name in actual_columns.items():
                    # Get the node attribute name from feature descriptions
                    node_attr = dssp_features[feature].split('(')[1].split(')')[0] if '(' in dssp_features[feature] else feature

                    # Add the feature to the node
                    if col_name in matching_rows.columns:
                        data[node_attr] = matching_rows[col_name].values[0]
                        features_added.add(node_attr)

                updated_nodes += 1

            except (ValueError, KeyError) as e:
                # Skip this node if there are issues
                continue

        print(f"Updated {updated_nodes} out of {len(graph.nodes)} nodes")
        print(f"Added the following features to nodes: {sorted(list(features_added))}")

        return graph
