import pytest
import numpy as np
import networkx as nx
import torch
from unittest.mock import patch, MagicMock
from torch_geometric.data import Data


from ProteinFeatureProcessor import ProteinFeatureProcessor

@pytest.fixture
def processor():
    """Fixture to create a ProteinFeatureProcessor instance"""
    return ProteinFeatureProcessor()

@pytest.fixture
def sample_graph():
    """Fixture to create a sample graph for testing"""
    graph = nx.Graph()

    # Add nodes with different features
    graph.add_node("A:1",
                   residue_name="ALA",
                   ss="H",
                   coords=np.array([1.0, 2.0, 3.0]),
                   b_factor=10.0,
                   meiler=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]))

    graph.add_node("A:2",
                   residue_name="CYS",
                   ss="E",
                   coords=np.array([4.0, 5.0, 6.0]),
                   b_factor=20.0,
                   meiler=np.array([0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]))

    # Add a node with missing features
    graph.add_node("A:3",
                   residue_name="UNK")

    # Add a node with non-standard residue
    graph.add_node("A:4",
                   residue_name="XYZ",  # Non-standard residue
                   ss="?",
                   coords=np.array([7.0, 8.0, 9.0]),
                   b_factor=30.0,
                   meiler=np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]))

    # Add edges with different features
    graph.add_edge("A:1", "A:2", distance=3.8, kind={"peptide_bond"})
    graph.add_edge("A:2", "A:3", distance=7.2, kind={"contact"})
    graph.add_edge("A:3", "A:4", distance=5.5)  # No kind specified

    return graph

def test_init(processor):
    """Test initialization of ProteinFeatureProcessor"""
    # Check amino acid vocabulary and mapping
    assert len(processor.amino_acids) == 21  # 20 standard + UNK
    assert len(processor.aa_to_idx) == 21
    assert processor.aa_to_idx["ALA"] == 0
    assert processor.aa_to_idx["UNK"] == 20

    # Check secondary structure vocabulary and mapping
    assert len(processor.ss_types) == 7
    assert len(processor.ss_to_idx) == 7
    assert processor.ss_to_idx["H"] == 0
    assert processor.ss_to_idx["?"] == 6

def test_process_graph_to_pyg(processor, sample_graph, monkeypatch):
    """Test processing a graph to PyTorch Geometric format"""
    # Mock the extract methods to isolate the test
    mock_extract_nodes = MagicMock(return_value=(
        np.array([[1.0, 2.0], [3.0, 4.0]]),  # Node features
        {"A:1": 0, "A:2": 1}  # Node mapping
    ))

    mock_extract_edges = MagicMock(return_value=(
        np.array([[0, 1], [1, 0]]).T,  # Edge indices
        np.array([[0.5], [0.6]])  # Edge features
    ))

    monkeypatch.setattr(processor, 'extract_node_features', mock_extract_nodes)
    monkeypatch.setattr(processor, 'extract_edge_features', mock_extract_edges)

    # Process the graph
    data = processor.process_graph_to_pyg(sample_graph)

    # Check the result
    assert isinstance(data, Data)
    assert data.x.shape == (2, 2)  # 2 nodes, 2 features
    assert data.edge_index.shape == (2, 2)  # 2 edges
    assert data.edge_attr.shape == (2, 1)  # 2 edges, 1 feature

    # Check the method calls
    mock_extract_nodes.assert_called_once_with(sample_graph)
    mock_extract_edges.assert_called_once()

def test_extract_node_features(processor, sample_graph):
    """Test extraction of node features"""
    node_features, node_mapping = processor.extract_node_features(sample_graph)

    # Check dimensions
    num_nodes = len(sample_graph.nodes)
    feature_dim = len(processor.amino_acids) + len(processor.ss_types) + 3 + 1 + 7
    assert node_features.shape == (num_nodes, feature_dim)

    # Check node mapping
    assert len(node_mapping) == num_nodes
    for node in sample_graph.nodes:
        assert node in node_mapping

    # Check that one-hot encoding works for amino acids
    for i, node in enumerate(sample_graph.nodes):
        residue_name = sample_graph.nodes[node].get('residue_name', 'UNK')
        if residue_name in processor.aa_to_idx:
            aa_idx = processor.aa_to_idx[residue_name]
        else:
            aa_idx = processor.aa_to_idx['UNK']
        assert node_features[i, aa_idx] == 1.0

    # Check that missing values are handled properly
    missing_node_idx = node_mapping["A:3"]
    ss_start_idx = len(processor.amino_acids)
    unk_ss_idx = processor.ss_to_idx['?']
    assert node_features[missing_node_idx, ss_start_idx + unk_ss_idx] == 1.0

    # Check that non-standard residues are mapped to UNK
    nonstandard_node_idx = node_mapping["A:4"]
    unk_aa_idx = processor.aa_to_idx['UNK']
    assert node_features[nonstandard_node_idx, unk_aa_idx] == 1.0

def test_extract_edge_features(processor, sample_graph):
    """Test extraction of edge features"""
    # First get the node mapping
    _, node_mapping = processor.extract_node_features(sample_graph)

    # Then extract edge features
    edge_indices, edge_features = processor.extract_edge_features(sample_graph, node_mapping)

    # Check dimensions
    num_edges = len(sample_graph.edges)
    assert edge_indices.shape == (2, num_edges)  # Shape should be (2, num_edges)
    assert edge_features.shape == (num_edges, 2)  # 2 features: distance and edge_type

    # Check edge indices
    for i, (u, v) in enumerate(sample_graph.edges):
        u_idx, v_idx = node_mapping[u], node_mapping[v]
        assert (edge_indices[:, i] == [u_idx, v_idx]).all() or (edge_indices[:, i] == [v_idx, u_idx]).all()

    # Check peptide bond feature
    for i, (u, v, data) in enumerate(sample_graph.edges(data=True)):
        edge_type_value = edge_features[i, 1]
        if 'kind' in data and 'peptide_bond' in data['kind']:
            assert edge_type_value == 1.0  # Peptide bond flag
        else:
            assert edge_type_value == 0.0  # Not a peptide bond

@pytest.mark.parametrize("features,coord_indices,b_factor_idx,expected_shapes", [
    # Simple case with coordinates and b-factor
    (
            np.array([
                [1.0, 2.0, 3.0, 10.0, 0.1, 0.2],
                [4.0, 5.0, 6.0, 20.0, 0.3, 0.4],
                [7.0, 8.0, 9.0, 30.0, 0.5, 0.6]
            ]),
            [0, 1, 2],
            3,
            (3, 6)
    ),
    # Case with zero standard deviation in coordinates
    (
            np.array([
                [1.0, 2.0, 3.0, 10.0],
                [1.0, 2.0, 3.0, 20.0],
                [1.0, 2.0, 3.0, 30.0]
            ]),
            [0, 1, 2],
            3,
            (3, 4)
    ),
])
def test_normalize_features(processor, features, coord_indices, b_factor_idx, expected_shapes):
    """Test normalization of features with different inputs"""
    normalized = processor.normalize_features(features, coord_indices, b_factor_idx)

    # Check shape
    assert normalized.shape == expected_shapes

    # Check that the original array wasn't modified
    assert not np.array_equal(features, normalized)

    # Check coordinate normalization
    coords = features[:, coord_indices]
    coords_mean = np.mean(coords, axis=0)
    coords_std = np.std(coords, axis=0)
    coords_std[coords_std == 0] = 1.0  # Avoid division by zero
    expected_coords = (coords - coords_mean) / coords_std
    np.testing.assert_array_almost_equal(normalized[:, coord_indices], expected_coords)

    # Check B-factor normalization
    b_factors = features[:, b_factor_idx]
    b_mean = np.mean(b_factors)
    b_std = np.std(b_factors)
    if b_std > 0:
        expected_b = (b_factors - b_mean) / b_std
        np.testing.assert_array_almost_equal(normalized[:, b_factor_idx], expected_b)

    # Check that other features weren't modified
    if b_factor_idx + 1 < features.shape[1]:
        np.testing.assert_array_equal(normalized[:, b_factor_idx+1:], features[:, b_factor_idx+1:])

# Add this block only if running directly
if __name__ == "__main__":
    pytest.main(["-xvs", __file__])