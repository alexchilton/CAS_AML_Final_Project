from graphein.protein.config import ProteinGraphConfig, DSSPConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.features.nodes.dssp import add_dssp_df
from graphein.protein.utils import download_pdb
import graphein.protein as gp
from functools import partial
import networkx as nx
import numpy as np
import gc
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



    def calculate_torsions_memory_efficient(self, chain_atoms, positions, distance_cutoff=2.0, backbone_only=True):
        """
        Calculate torsions with memory efficiency.

        Parameters:
        -----------
        chain_atoms : list
            List of atom data
        positions : ndarray
            Array of atom positions
        distance_cutoff : float
            Maximum distance to consider for connections
        backbone_only : bool
            If True, only calculate torsions for backbone atoms

        Returns:
        --------
        dict
            Dictionary of torsion angles
        """
        if backbone_only:
            return self.calculate_backbone_torsions(chain_atoms, positions)
        else:
            return self.calculate_all_torsions(chain_atoms, positions, distance_cutoff)


    def calculate_backbone_torsions(self, chain_atoms, positions):
        """
        Calculate only backbone torsion angles (phi, psi, omega) for a protein chain.

        Parameters:
        -----------
        chain_atoms : list
            List of atom data
        positions : ndarray
            Array of atom positions

        Returns:
        --------
        dict
            Dictionary of backbone torsion angles
        """
        torsions = {}

        # Group atoms by residue
        residue_groups = {}

        # First pass: organize atoms by residue
        for i, atom in enumerate(chain_atoms):
            chain_id, res_id, res_name, atom_name, element, position = atom

            if atom_name in ['N', 'CA', 'C', 'O']:
                if res_id not in residue_groups:
                    residue_groups[res_id] = {}
                residue_groups[res_id][atom_name] = i

        # Sort residues by ID
        sorted_res_ids = sorted(residue_groups.keys())

        # Calculate phi, psi, omega angles between consecutive residues
        for i in range(1, len(sorted_res_ids)):
            prev_res = sorted_res_ids[i-1]
            curr_res = sorted_res_ids[i]

            # Skip if residues aren't sequential
            if curr_res - prev_res != 1:
                continue

            prev_atoms = residue_groups.get(prev_res, {})
            curr_atoms = residue_groups.get(curr_res, {})

            # Calculate phi (C-N-CA-C)
            phi_torsion = self.calculate_phi_torsion(prev_atoms, curr_atoms, positions)
            if phi_torsion:
                torsions[phi_torsion[0]] = phi_torsion[1]

            # Calculate psi if next residue exists
            if i < len(sorted_res_ids) - 1:
                next_res = sorted_res_ids[i+1]

                # Skip if residues aren't sequential
                if next_res - curr_res != 1:
                    continue

                next_atoms = residue_groups.get(next_res, {})

                psi_torsion = self.calculate_psi_torsion(curr_atoms, next_atoms, positions)
                if psi_torsion:
                    torsions[psi_torsion[0]] = psi_torsion[1]

            # Calculate omega
            omega_torsion = self.calculate_omega_torsion(prev_atoms, curr_atoms, positions)
            if omega_torsion:
                torsions[omega_torsion[0]] = omega_torsion[1]

        print(f"    Calculated {len(torsions)} backbone torsions")
        return torsions


    def calculate_phi_torsion(self, prev_atoms, curr_atoms, positions):
        """
        Calculate phi torsion angle (C-N-CA-C).

        Parameters:
        -----------
        prev_atoms : dict
            Dictionary of previous residue atoms by name->index
        curr_atoms : dict
            Dictionary of current residue atoms by name->index
        positions : ndarray
            Array of atom positions

        Returns:
        --------
        tuple or None
            ((i, j, k, l), torsion_value) or None if calculation fails
        """
        if ('C' in prev_atoms and 'N' in curr_atoms and
                'CA' in curr_atoms and 'C' in curr_atoms):
            i = prev_atoms['C']
            j = curr_atoms['N']
            k = curr_atoms['CA']
            l = curr_atoms['C']

            torsion = self.calculate_single_torsion(
                positions[i], positions[j], positions[k], positions[l])
            if torsion is not None:
                return ((i, j, k, l), torsion)
        return None


    def calculate_psi_torsion(self,curr_atoms, next_atoms, positions):
        """
        Calculate psi torsion angle (N-CA-C-N).

        Parameters:
        -----------
        curr_atoms : dict
            Dictionary of current residue atoms by name->index
        next_atoms : dict
            Dictionary of next residue atoms by name->index
        positions : ndarray
            Array of atom positions

        Returns:
        --------
        tuple or None
            ((i, j, k, l), torsion_value) or None if calculation fails
        """
        if ('N' in curr_atoms and 'CA' in curr_atoms and
                'C' in curr_atoms and 'N' in next_atoms):
            i = curr_atoms['N']
            j = curr_atoms['CA']
            k = curr_atoms['C']
            l = next_atoms['N']

            torsion = self.calculate_single_torsion(
                positions[i], positions[j], positions[k], positions[l])
            if torsion is not None:
                return ((i, j, k, l), torsion)
        return None


    def calculate_omega_torsion(self,prev_atoms, curr_atoms, positions):
        """
        Calculate omega torsion angle (CA-C-N-CA).

        Parameters:
        -----------
        prev_atoms : dict
            Dictionary of previous residue atoms by name->index
        curr_atoms : dict
            Dictionary of current residue atoms by name->index
        positions : ndarray
            Array of atom positions

        Returns:
        --------
        tuple or None
            ((i, j, k, l), torsion_value) or None if calculation fails
        """
        if ('CA' in prev_atoms and 'C' in prev_atoms and
                'N' in curr_atoms and 'CA' in curr_atoms):
            i = prev_atoms['CA']
            j = prev_atoms['C']
            k = curr_atoms['N']
            l = curr_atoms['CA']

            torsion = self.calculate_single_torsion(
                positions[i], positions[j], positions[k], positions[l])
            if torsion is not None:
                return ((i, j, k, l), torsion)
        return None


    def calculate_all_torsions(self, chain_atoms, positions, distance_cutoff=2.0):
        """
        Calculate all possible torsion angles for a protein chain.

        Parameters:
        -----------
        chain_atoms : list
            List of atom data
        positions : ndarray
            Array of atom positions
        distance_cutoff : float
            Maximum distance to consider for connections

        Returns:
        --------
        dict
            Dictionary of torsion angles
        """
        torsions = {}
        n_atoms = len(chain_atoms)

        # First identify connected pairs (similar to angle calculation)
        connected_pairs = self.identify_connected_pairs(chain_atoms, positions, distance_cutoff)

        # Find valid quadruplets (i-j-k-l where each adjacent pair is connected)
        jk_pairs = self.identify_connected_jk_pairs(connected_pairs, n_atoms)

        # Process j-k pairs in batches
        torsions = self.process_torsion_batches(jk_pairs, connected_pairs, positions, n_atoms)

        print(f"    Calculated {len(torsions)} torsions")
        return torsions


    def identify_connected_pairs(slef, chain_atoms, positions, distance_cutoff):
        """
        Identify pairs of atoms that are connected.

        Parameters:
        -----------
        chain_atoms : list
            List of atom data
        positions : ndarray
            Array of atom positions
        distance_cutoff : float
            Maximum distance to consider for connections

        Returns:
        --------
        set
            Set of connected atom pairs (i,j)
        """
        connected_pairs = set()
        n_atoms = len(chain_atoms)

        for i in range(n_atoms):
            res_i = chain_atoms[i][1]
            pos_i = positions[i]

            for j in range(i+1, n_atoms):
                res_j = chain_atoms[j][1]
                pos_j = positions[j]

                # Check if connected by residue
                residue_connected = (res_i == res_j or abs(res_i - res_j) == 1)

                # Check distance
                distance = np.linalg.norm(pos_j - pos_i)

                if residue_connected or distance <= distance_cutoff:
                    connected_pairs.add((i, j))
                    connected_pairs.add((j, i))  # Undirected

        return connected_pairs


    def identify_connected_jk_pairs(self, connected_pairs, n_atoms):
        """
        Identify all connected j-k pairs for torsion calculation.

        Parameters:
        -----------
        connected_pairs : set
            Set of connected atom pairs
        n_atoms : int
            Total number of atoms

        Returns:
        --------
        list
            List of j-k pairs where j < k
        """
        jk_pairs = []
        for j in range(n_atoms):
            for k in range(j+1, n_atoms):
                if (j, k) in connected_pairs:
                    jk_pairs.append((j, k))

        return jk_pairs


    def process_torsion_batches(self, jk_pairs, connected_pairs, positions, n_atoms, batch_size=500):
        """
        Process torsion calculations in batches for memory efficiency.

        Parameters:
        -----------
        jk_pairs : list
            List of j-k pairs
        connected_pairs : set
            Set of connected atom pairs
        positions : ndarray
            Array of atom positions
        n_atoms : int
            Total number of atoms
        batch_size : int
            Size of batches for processing

        Returns:
        --------
        dict
            Dictionary of torsion angles
        """
        torsions = {}
        quadruplet_count = 0

        # Process j-k pairs in batches
        for jk_start in range(0, len(jk_pairs), batch_size):
            jk_end = min(jk_start + batch_size, len(jk_pairs))

            for jk_idx in range(jk_start, jk_end):
                j, k = jk_pairs[jk_idx]

                # Find i's connected to j
                i_candidates = self.find_connected_atoms(j, k, connected_pairs, n_atoms)

                # Find l's connected to k
                l_candidates = self.find_connected_atoms(k, j, connected_pairs, n_atoms)

                # Process i-l candidates in smaller batches
                quadruplet_count = self.process_ijkl_candidates(
                    i_candidates, l_candidates, j, k, positions, torsions,
                    quadruplet_count, batch_size
                )

            # Force garbage collection
            gc.collect()

        return torsions


    def find_connected_atoms(self, atom_idx, exclude_idx, connected_pairs, n_atoms):
        """
        Find atoms connected to a specific atom.

        Parameters:
        -----------
        atom_idx : int
            Index of the atom to find connections for
        exclude_idx : int
            Index of atom to exclude
        connected_pairs : set
            Set of connected atom pairs
        n_atoms : int
            Total number of atoms

        Returns:
        --------
        list
            List of atoms connected to atom_idx
        """
        candidates = []
        for idx in range(n_atoms):
            if idx != atom_idx and idx != exclude_idx and (idx, atom_idx) in connected_pairs:
                candidates.append(idx)
        return candidates


    def process_ijkl_candidates(self, i_candidates, l_candidates, j, k, positions, torsions,
                                quadruplet_count, batch_size):
        """
        Process i-j-k-l candidates for torsion calculation.

        Parameters:
        -----------
        i_candidates : list
            List of atoms connected to j
        l_candidates : list
            List of atoms connected to k
        j : int
            Index of atom j
        k : int
            Index of atom k
        positions : ndarray
            Array of atom positions
        torsions : dict
            Dictionary to store calculated torsions
        quadruplet_count : int
            Running count of calculated torsions
        batch_size : int
            Size of batches for processing

        Returns:
        --------
        int
            Updated count of calculated torsions
        """
        for i_start in range(0, len(i_candidates), batch_size):
            i_end = min(i_start + batch_size, len(i_candidates))

            for l_start in range(0, len(l_candidates), batch_size):
                l_end = min(l_start + batch_size, len(l_candidates))

                for i_idx in range(i_start, i_end):
                    i = i_candidates[i_idx]

                    for l_idx in range(l_start, l_end):
                        l = l_candidates[l_idx]

                        # Skip if i and l are the same atom
                        if i == l:
                            continue

                        # Calculate torsion
                        torsion = self.calculate_single_torsion(
                            positions[i], positions[j], positions[k], positions[l])
                        if torsion is not None:
                            torsions[(i, j, k, l)] = torsion
                            quadruplet_count += 1

                # Force garbage collection
                gc.collect()

        return quadruplet_count


    def calculate_single_torsion(self, pos_i, pos_j, pos_k, pos_l):
        """
        Calculate a single torsion angle from 4 positions.

        Parameters:
        -----------
        pos_i, pos_j, pos_k, pos_l : ndarray
            3D coordinates of the four atoms

        Returns:
        --------
        float or None
            Torsion angle in degrees or None if calculation fails
        """
        # Define vectors
        v1 = pos_j - pos_i
        v2 = pos_k - pos_j
        v3 = pos_l - pos_k

        # Skip if any vector is too short
        if (np.linalg.norm(v1) < 1e-6 or
                np.linalg.norm(v2) < 1e-6 or
                np.linalg.norm(v3) < 1e-6):
            return None

        # Cross products
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        # Skip if either normal is too short
        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)

        if n1_norm < 1e-6 or n2_norm < 1e-6:
            return None

        # Normalize
        n1 = n1 / n1_norm
        n2 = n2 / n2_norm

        # Compute torsion angle
        v2_norm = np.linalg.norm(v2)
        torsion = np.arctan2(
            np.dot(np.cross(n1, n2), v2 / v2_norm),
            np.dot(n1, n2)
        )

        return np.degrees(torsion)
