from Bio.PDB.PDBParser import PDBParser
import numpy as np
import pandas as pd
import torch
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
from Bio.PDB.PDBParser import PDBParser
from sklearn.preprocessing import LabelEncoder
from torch_geometric.data import Data
import networkx as nx
import matplotlib.pyplot as plt

### would be better if everything was being read in as floats...!!!


class ProteinAnalyzer:
    def __init__(self, pdb_file, aa_info_file):
        self.pdb_file = pdb_file
        self.aa_info_file = aa_info_file
        self.parser = PDBParser(PERMISSIVE=1)
        self.structure = self.parser.get_structure('protein_structure', self.pdb_file)
        self.aa_info_dict = self.get_aa_info_dict()
        self.aa_info_dict_short = self.get_aa_info_dict_short()
        self.c_alpha_df = self.extract_c_alpha_info()

    def get_aa_info_dict(self):
        '''Read the amino acid info file and return a dictionary with the amino acid properties.'''
        aa_info = pd.read_csv(self.aa_info_file)
        return aa_info.set_index('Abbrev.').to_dict(orient='index')

    def get_aa_info_dict_short(self):
        '''Read the amino acid info file and return a dictionary with the amino acid properties.'''
        aa_info = pd.read_csv(self.aa_info_file)
        return aa_info.set_index('Short').to_dict(orient='index')

    def get_residue_info(self, res_name):
        '''Return the short name of the amino acid given the residue name.'''
        if res_name in self.aa_info_dict:
            return self.aa_info_dict[res_name]['Short']
        else:
            print(f"Residue {res_name} not in amino acid info dictionary.")
            return None

    def process_residue(self, residue, res_name, c_alpha_matrix):
        '''Process a residue and extract the C-alpha atom coordinates,.'''
        if 'CA' in residue:
            ca_atom = residue['CA']
            x, y, z = ca_atom.get_coord()
            c_alpha_matrix.append([
                x, y, z,
                self.aa_info_dict[res_name]['Short'],
                self.aa_info_dict[res_name]['Avg. mass (Da)']
            ])
        else:
            print(f"No C-alpha atom found for residue {res_name} in chain {residue.parent.id}")

    def create_c_alpha_dataframe(self, c_alpha_matrix):
        '''Create a DataFrame from the C-alpha matrix.'''
        c_alpha_df = pd.DataFrame(c_alpha_matrix, columns=["X", "Y", "Z", "AA", "Mass"])
        if c_alpha_df.empty:
            print("C-alpha DataFrame is empty. Check residue names or PDB file.")
        #else:
        #print(c_alpha_df)
        return c_alpha_df

    def extract_c_alpha_info(self):
        '''Extract the C-alpha atom coordinates and amino acid properties from the PDB file.'''
        c_alpha_matrix = []
        res_list = Selection.unfold_entities(self.structure, "R")
        for residue in res_list:
            res_name = residue.get_resname().title()
            self.process_residue(residue, res_name, c_alpha_matrix)
        return self.create_c_alpha_dataframe(c_alpha_matrix)

    def calculate_pairwise_distances(self, coordinates):
        '''Calculate the pairwise distances between C-alpha atoms.'''
        num_residues = len(coordinates)
        distance_matrix = np.zeros((num_residues, num_residues))
        for i in range(num_residues):
            for j in range(i, num_residues):
                dist = np.linalg.norm(coordinates[i] - coordinates[j])
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist
        return distance_matrix

    def encode_amino_acid_properties(self, aa_letters):
        '''Encode the amino acid properties for the given amino acid letters.'''
        encoded_features = []
        for aa in aa_letters:
            aa_str = aa[0]  # Convert numpy array to string
            if aa_str in self.aa_info_dict_short:
                avg_mass = self.aa_info_dict_short[aa_str]['Avg. mass (Da)']
                encoded_features.append([aa_str, avg_mass])
            else:
                encoded_features.append([aa_str, 0])
        return np.array(encoded_features)

    def prepare_autoencoder_input(self):
        '''Prepare the input DataFrame for the autoencoder.'''
        coords = self.c_alpha_df[['X', 'Y', 'Z']].values
        aa_letters = self.c_alpha_df['AA'].values.reshape(-1, 1)  # Reshape to add as a column
        encoded_features = self.encode_amino_acid_properties(aa_letters)
        neighborhood_info = self.calculate_neighborhood_info()
        autoencoder_input = np.hstack([coords, encoded_features, np.array(neighborhood_info)])
        columns = ['X', 'Y', 'Z', 'AA','Avg_Mass', 'Avg_Neighbor_Dist', 'Max_Neighbor_Dist', 'Neighbor_Count']
        return pd.DataFrame(autoencoder_input, columns=columns)


    def calculate_neighborhood_info(self, neighborhood_radius=5.0):
        '''Calculate the neighborhood information for each C-alpha atom.'''
        atom_list = [atom for atom in self.structure.get_atoms() if atom.name == 'CA']
        neighbor_search = NeighborSearch(atom_list)
        neighborhood_info = []
        for res in self.structure.get_residues():
            if 'CA' in res:
                center = res['CA'].get_coord()
                neighbors = neighbor_search.search(center, neighborhood_radius)
                distances = [np.linalg.norm(center - neighbor.get_coord()) for neighbor in neighbors if neighbor != res['CA']]
                avg_distance = np.mean(distances) if distances else 0
                max_distance = np.max(distances) if distances else 0
                count_neighbors = len(distances)
                neighborhood_info.append([avg_distance, max_distance, count_neighbors])
        return neighborhood_info




    def pad_dataframe(self, df, target_shape):
        '''Pad the DataFrame to the target shape by adding rows of zeros.
        NOT USED currently - i pad the graph'''
        current_shape = df.shape
        if current_shape[0] >= target_shape[0] and current_shape[1] == target_shape[1]:
            return df
        rows_to_add = target_shape[0] - current_shape[0]
        if rows_to_add > 0:
            padding_df = pd.DataFrame(np.zeros((rows_to_add, target_shape[1])), columns=df.columns)
            padded_df = pd.concat([df, padding_df], ignore_index=True)
        else:
            padded_df = df
        return padded_df

    def create_graph(df, distance_threshold=5.0):
        '''Create a graph Data object from the input DataFrame.'''
        # Label encode the 'AA' column
        le = LabelEncoder()
        df['AA_encoded'] = le.fit_transform(df['AA'])

        # Ensure all columns are numeric
        df[['X', 'Y', 'Z', 'Avg_Mass', 'Avg_Neighbor_Dist', 'Max_Neighbor_Dist', 'Neighbor_Count']] = df[['X', 'Y', 'Z', 'Avg_Mass', 'Avg_Neighbor_Dist', 'Max_Neighbor_Dist', 'Neighbor_Count']].apply(pd.to_numeric, errors='coerce')

        # Extract node features (XYZ coordinates, encoded AA, and other features)
        node_features = torch.tensor(df[['X', 'Y', 'Z', 'AA_encoded', 'Avg_Mass', 'Avg_Neighbor_Dist', 'Max_Neighbor_Dist', 'Neighbor_Count']].values, dtype=torch.float)

        # Calculate pairwise distances
        coordinates = df[['X', 'Y', 'Z']].values
        num_nodes = coordinates.shape[0]
        distance_matrix = np.zeros((num_nodes, num_nodes))
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):
                dist = np.linalg.norm(coordinates[i] - coordinates[j])
                distance_matrix[i, j] = dist
                distance_matrix[j, i] = dist

        # Define edges based on the distance threshold
        edge_index = []
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):
                if distance_matrix[i, j] <= distance_threshold:
                    edge_index.append([i, j])
                    edge_index.append([j, i])

        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()

        # Create the Data object
        data = Data(x=node_features, edge_index=edge_index)
        data.aa_labels = df['AA'].values  # Store the amino acid labels
        return data, le

    def print_graph_metrics(graph):
        '''Print the number of nodes, edges, and the shape of the node features.
        not sure if used anywhere...duplicate of the one in GraphCreatorOneHotEncodedVariant'''
        print(f"Number of nodes: {graph.num_nodes}")
        print(f"Number of edges: {graph.num_edges}")
        print(f"Node features shape: {graph.x.shape}")

    def draw_graph(graph, le):
        '''Draw the graph using NetworkX.'''
        G = nx.Graph()
        edge_index = graph.edge_index.numpy()
        for i in range(edge_index.shape[1]):
            G.add_edge(edge_index[0, i], edge_index[1, i])

        pos = nx.spring_layout(G)
        aa_labels = graph.aa_labels
        node_colors = [plt.cm.tab20(le.transform([aa])[0] % 20) for aa in aa_labels]  # Color code nodes

        nx.draw(G, pos, with_labels=True, labels={i: aa_labels[i] for i in range(len(aa_labels))}, node_color=node_colors, node_size=50, font_size=8)
        plt.show()

# pdb_file = '/Users/alexchilton/Downloads/archive/train/AF-D0ZA02-F1-model_v4.pdb'
# aa_info_file = 'aa_mass_letter.csv'
# analyzer = ProteinAnalyzer(pdb_file, aa_info_file)
#
# # Generate the autoencoder input DataFrame
# autoencoder_input_df = analyzer.prepare_autoencoder_input()
#
# pd.set_option('display.max_columns', None)
#
# print(autoencoder_input_df.shape)
# print(autoencoder_input_df.head())