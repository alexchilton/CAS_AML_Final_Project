import pytest
import os
import networkx as nx
import torch
import pickle
from unittest.mock import patch, mock_open, MagicMock
from torch_geometric.data import Data

# Import the GraphBuilder class
from GraphBuilder import GraphBuilder

@pytest.fixture
def mock_base_dir():
    return "mock_base_dir"

@pytest.fixture
def mock_nanobody_dirs():
    return ["nb1", "nb2", "nb3"]

@pytest.fixture
def graph_builder(mock_base_dir, monkeypatch):
    # Patch os functions
    monkeypatch.setattr(os.path, 'exists', lambda x: True)
    monkeypatch.setattr(os, 'listdir', lambda x: ["nb1", "nb2", "nb3"])
    monkeypatch.setattr(os.path, 'isdir', lambda x: True)

    # Create GraphBuilder instance
    return GraphBuilder(mock_base_dir)

def test_init(graph_builder, mock_base_dir, mock_nanobody_dirs):
    """Test initialization of GraphBuilder"""
    assert graph_builder.base_dir == mock_base_dir
    assert graph_builder.nanobody_dirs == mock_nanobody_dirs

def test_find_nanobody_dirs(graph_builder, mock_nanobody_dirs, monkeypatch):
    """Test finding nanobody directories"""
    assert graph_builder._find_nanobody_dirs() == mock_nanobody_dirs

    # Test when base directory doesn't exist
    monkeypatch.setattr(os.path, 'exists', lambda x: False)
    with pytest.raises(FileNotFoundError):
        graph_builder._find_nanobody_dirs()

def test_get_nanobody_ids(graph_builder, mock_nanobody_dirs):
    """Test getting nanobody IDs"""
    assert graph_builder.get_nanobody_ids() == mock_nanobody_dirs

@pytest.mark.parametrize("nanobody_id,path_exists,expected_result,expected_exception", [
    ("nb1", True, "mock_graph", None),  # Success case
    ("nonexistent", True, None, ValueError),  # Nanobody not found
    ("nb1", False, None, FileNotFoundError),  # File not found
])
def test_load_graph(graph_builder, nanobody_id, path_exists, expected_result, expected_exception, monkeypatch):
    """Test loading a graph for a specific nanobody with parametrized cases"""
    # Mock pickle.load
    mock_graph = nx.Graph()

    # Set up monkeypatches
    monkeypatch.setattr(os.path, 'exists', lambda x: path_exists)
    monkeypatch.setattr(pickle, 'load', lambda f: mock_graph)
    monkeypatch.setattr(builtins, 'open', mock_open())

    if expected_exception:
        with pytest.raises(expected_exception):
            graph_builder.load_graph(nanobody_id)
    else:
        result = graph_builder.load_graph(nanobody_id)
        assert result == mock_graph

@pytest.mark.parametrize("nanobody_id,path_exists,expected_result,expected_exception", [
    ("nb1", True, {"key": "value"}, None),  # Success case
    ("nonexistent", True, None, ValueError),  # Nanobody not found
    ("nb1", False, None, FileNotFoundError),  # File not found
])
def test_load_data_dict(graph_builder, nanobody_id, path_exists, expected_result, expected_exception, monkeypatch):
    """Test loading data dictionary for a specific nanobody with parametrized cases"""
    # Mock data
    mock_data = {"key": "value"}

    # Set up monkeypatches
    monkeypatch.setattr(os.path, 'exists', lambda x: path_exists)
    monkeypatch.setattr(pickle, 'load', lambda f: mock_data)
    monkeypatch.setattr(builtins, 'open', mock_open())

    if expected_exception:
        with pytest.raises(expected_exception):
            graph_builder.load_data_dict(nanobody_id)
    else:
        result = graph_builder.load_data_dict(nanobody_id)
        assert result == mock_data

def test_list_chain_files(graph_builder, monkeypatch):
    """Test listing chain files for a nanobody"""
    # Mock directory listing
    mock_files = [
        "nb1_graph.pkl",
        "nb1_data.pkl",
        "nb1_A_atoms.pkl",
        "nb1_A_backbone.pkl",
        "nb1_B_atoms.pkl",
        "random_file.txt"
    ]
    monkeypatch.setattr(os, 'listdir', lambda x: mock_files)

    chain_files = graph_builder.list_chain_files("nb1")
    expected = {
        "A": ["nb1_A_atoms.pkl", "nb1_A_backbone.pkl"],
        "B": ["nb1_B_atoms.pkl"]
    }
    assert chain_files == expected

    # Test nanobody ID not found
    with pytest.raises(ValueError):
        graph_builder.list_chain_files("nonexistent")

@pytest.mark.parametrize("nanobody_id,chain_id,file_type,path_exists,expected_result,expected_exception", [
    ("nb1", "A", "atoms", True, {"chain_data": "value"}, None),  # Success case
    ("nonexistent", "A", "atoms", True, None, ValueError),  # Nanobody not found
    ("nb1", "A", "atoms", False, None, FileNotFoundError),  # File not found
])
def test_load_chain_file(graph_builder, nanobody_id, chain_id, file_type, path_exists,
                         expected_result, expected_exception, monkeypatch):
    """Test loading a specific chain file with parametrized cases"""
    # Mock data
    mock_data = {"chain_data": "value"}

    # Set up monkeypatches
    monkeypatch.setattr(os.path, 'exists', lambda x: path_exists)
    monkeypatch.setattr(pickle, 'load', lambda f: mock_data)
    monkeypatch.setattr(builtins, 'open', mock_open())

    if expected_exception:
        with pytest.raises(expected_exception):
            graph_builder.load_chain_file(nanobody_id, chain_id, file_type)
    else:
        result = graph_builder.load_chain_file(nanobody_id, chain_id, file_type)
        assert result == mock_data

def test_get_available_file_types(graph_builder, monkeypatch):
    """Test getting available file types for a chain"""
    # Mock list_chain_files
    mock_chain_files = {
        "A": ["nb1_A_atoms.pkl", "nb1_A_backbone.pkl"],
        "B": ["nb1_B_atoms.pkl"]
    }
    monkeypatch.setattr(graph_builder, 'list_chain_files', lambda x: mock_chain_files)

    file_types = graph_builder.get_available_file_types("nb1", "A")
    assert sorted(file_types) == sorted(["atoms", "backbone"])

    # Test chain ID not found
    with pytest.raises(ValueError):
        graph_builder.get_available_file_types("nb1", "nonexistent")

def test_print_graph_summary(capsys):
    """Test printing graph summary"""
    # Create a simple graph
    mock_graph = nx.Graph()
    mock_graph.add_node(1, attr="value")
    mock_graph.add_edge(1, 2, weight=0.5)

    # Create GraphBuilder instance
    builder = GraphBuilder("mock_dir")

    # Call print_graph_summary
    builder.print_graph_summary(mock_graph)

    # Capture output
    captured = capsys.readouterr()

    # Verify output contains expected information
    assert "Graph Summary" in captured.out
    assert "Number of nodes: 2" in captured.out
    assert "Number of edges: 1" in captured.out

def test_convert_to_pytorch_geometric():
    """Test converting NetworkX graph to PyTorch Geometric Data"""
    # Create a mock graph
    mock_graph = nx.Graph()
    mock_graph.add_node("node1", coords=[1.0, 2.0, 3.0])
    mock_graph.add_node("node2", CA=[4.0, 5.0, 6.0])
    mock_graph.add_node("node3")  # No coords
    mock_graph.add_edge("node1", "node2", distance=2.5)
    mock_graph.add_edge("node2", "node3")  # No distance

    # Create GraphBuilder instance
    builder = GraphBuilder("mock_dir")

    # Convert to PyG
    data = builder.convert_to_pytorch_geometric(mock_graph)

    # Check the result
    assert isinstance(data, Data)
    assert data.x.shape == (3, 3)  # 3 nodes, 3 features (coords)
    assert data.edge_index.shape == (2, 2)  # 2 edges
    assert data.edge_attr.shape == (2, 1)  # 2 edges, 1 feature (distance)

import builtins  # Need to import to monkeypatch

# Add this block only if running directly
if __name__ == "__main__":
    pytest.main(["-xvs", __file__])