import pytest
import os
import numpy as np
import networkx as nx
import torch
from unittest.mock import patch, MagicMock, mock_open
from torch_geometric.data import Data

# Import the classes
from ProteinGraphDataset import ProteinGraphDataset
from GraphBuilder import GraphBuilder
from ProteinFeatureProcessor import ProteinFeatureProcessor

@pytest.fixture
def mock_graph_builder():
    """Create a mock GraphBuilder"""
    mock = MagicMock(spec=GraphBuilder)
    mock.get_nanobody_ids.return_value = ["nb1", "nb2", "nb3"]
    return mock

@pytest.fixture
def mock_feature_processor():
    """Create a mock ProteinFeatureProcessor"""
    return MagicMock(spec=ProteinFeatureProcessor)

@pytest.fixture
def dataset(monkeypatch, mock_graph_builder, mock_feature_processor):
    """Create a ProteinGraphDataset with mocked dependencies"""
    # Patch os functions
    monkeypatch.setattr(os.path, 'exists', lambda x: True)
    monkeypatch.setattr(os, 'makedirs', lambda x: None)

    # Patch class constructors
    monkeypatch.setattr('ProteinGraphDataset.GraphBuilder', lambda x: mock_graph_builder)
    monkeypatch.setattr('ProteinGraphDataset.ProteinFeatureProcessor', lambda: mock_feature_processor)

    # Create dataset
    return ProteinGraphDataset(root="mock_root_dir", use_cache=True)

def test_init(dataset, mock_graph_builder):
    """Test initialization of ProteinGraphDataset"""
    # Check basic properties
    assert dataset.root_dir == "mock_root_dir"
    assert dataset.use_cache is True
    assert dataset.selected_nanobodies is None

    # Check that nanobody_ids were retrieved from GraphBuilder
    mock_graph_builder.get_nanobody_ids.assert_called_once()
    assert dataset.nanobody_ids == ["nb1", "nb2", "nb3"]

def test_init_with_selected_nanobodies(monkeypatch, mock_graph_builder, mock_feature_processor):
    """Test initialization with selected nanobodies"""
    # Patch os functions
    monkeypatch.setattr(os.path, 'exists', lambda x: True)
    monkeypatch.setattr(os, 'makedirs', lambda x: None)

    # Patch class constructors
    monkeypatch.setattr('ProteinGraphDataset.GraphBuilder', lambda x: mock_graph_builder)
    monkeypatch.setattr('ProteinGraphDataset.ProteinFeatureProcessor', lambda: mock_feature_processor)

    # Create dataset with selected nanobodies
    selected = ["nb1", "nb3"]
    dataset = ProteinGraphDataset(
        root="mock_root_dir",
        use_cache=True,
        selected_nanobodies=selected
    )

    assert dataset.nanobody_ids == ["nb1", "nb3"]

def test_processed_dir(dataset):
    """Test processed_dir property"""
    expected_dir = os.path.join("mock_root_dir", 'processed')
    assert dataset.processed_dir == expected_dir

@pytest.mark.parametrize("all_ids,selected,expected", [
    (["nb1", "nb2", "nb3", "nb4"], None, ["nb1", "nb2", "nb3", "nb4"]),  # No filtering
    (["nb1", "nb2", "nb3", "nb4"], ["nb1", "nb3"], ["nb1", "nb3"]),  # With filtering
    (["nb1", "nb2", "nb3", "nb4"], ["nb5"], []),  # No matches
])
def test_filter_nanobodies(dataset, all_ids, selected, expected):
    """Test filtering of nanobody IDs with different inputs"""
    dataset.selected_nanobodies = selected
    filtered = dataset._filter_nanobodies(all_ids)
    assert filtered == expected

def test_raw_file_names(dataset):
    """Test raw_file_names property"""
    assert dataset.raw_file_names == []

def test_processed_file_names(dataset):
    """Test processed_file_names property"""
    expected = ["nb1.pt", "nb2.pt", "nb3.pt"]
    assert dataset.processed_file_names == expected

def test_len(dataset):
    """Test len method"""
    assert len(dataset) == 3

    # Test with empty dataset
    dataset.nanobody_ids = []
    assert len(dataset) == 0

@pytest.mark.parametrize("idx,nanobody_id", [
    (0, "nb1"),
    (1, "nb2"),
    (2, "nb3"),
])
def test_get_from_cache(dataset, idx, nanobody_id, monkeypatch):
    """Test get method with cached data"""
    # Set up the mock return for torch.load
    mock_data = Data(x=torch.tensor([1.0, 2.0]), edge_index=torch.tensor([[0], [1]]))

    # Configure patch for torch.load
    monkeypatch.setattr(torch, 'load', lambda x: mock_data)

    # Configure path exists to return True (cache exists)
    monkeypatch.setattr(os.path, 'exists', lambda x: True)

    # Get the data
    data = dataset.get(idx)

    # Check the result
    assert data == mock_data

@pytest.mark.parametrize("idx,nanobody_id", [
    (0, "nb1"),
    (1, "nb2"),
    (2, "nb3"),
])
def test_get_uncached(dataset, idx, nanobody_id, mock_graph_builder, mock_feature_processor, monkeypatch):
    """Test get method with uncached data"""
    # Configure path exists to return False (no cache)
    monkeypatch.setattr(os.path, 'exists', lambda x: False)

    # Set up mocks for processing chain
    mock_graph = nx.Graph()
    mock_data = Data(x=torch.tensor([3.0, 4.0]), edge_index=torch.tensor([[0], [1]]))

    mock_graph_builder.load_graph.return_value = mock_graph
    mock_feature_processor.process_graph_to_pyg.return_value = mock_data

    # Configure patch for torch.save
    monkeypatch.setattr(torch, 'save', lambda data, path: None)

    # Get the data
    with patch.object(dataset, '_process_nanobody', wraps=dataset._process_nanobody) as mock_process:
        data = dataset.get(idx)

        # Check that _process_nanobody was called
        mock_process.assert_called_once_with(nanobody_id)

    # Check the result
    assert data == mock_data

    # Check that the data was processed correctly
    mock_graph_builder.load_graph.assert_called_once_with(nanobody_id)
    mock_feature_processor.process_graph_to_pyg.assert_called_once_with(mock_graph)

def test_process_nanobody_success(dataset, mock_graph_builder, mock_feature_processor):
    """Test _process_nanobody method with successful processing"""
    # Set up successful processing
    mock_graph = nx.Graph()
    mock_data = Data(x=torch.tensor([5.0, 6.0]), edge_index=torch.tensor([[0], [1]]))

    mock_graph_builder.load_graph.return_value = mock_graph
    mock_feature_processor.process_graph_to_pyg.return_value = mock_data

    # Process a nanobody
    data = dataset._process_nanobody("nb1")

    # Check the result
    assert data == mock_data
    assert data.nanobody_id == "nb1"

    # Check that the processing chain was called correctly
    mock_graph_builder.load_graph.assert_called_once_with("nb1")
    mock_feature_processor.process_graph_to_pyg.assert_called_once_with(mock_graph)

def test_process_nanobody_error(dataset, mock_graph_builder, mock_feature_processor):
    """Test _process_nanobody method error handling"""
    # Set up error during processing
    mock_graph_builder.load_graph.side_effect = Exception("Test error")

    # Process a nanobody
    data = dataset._process_nanobody("nb1")

    # Check the result (should be an empty Data object)
    assert isinstance(data, Data)
    assert data.nanobody_id == "nb1"
    assert data.x.shape[0] == 0
    assert data.edge_index.shape[1] == 0

def test_process(dataset, monkeypatch):
    """Test process method"""
    # Mock _process_nanobody and torch.save
    process_mock = MagicMock(return_value=Data())
    monkeypatch.setattr(dataset, '_process_nanobody', process_mock)
    monkeypatch.setattr(torch, 'save', lambda data, path: None)

    # Mock os.path.exists for processed files check
    monkeypatch.setattr(os.path, 'exists', lambda path: False)

    # Call process
    dataset.process()

    # Check that _process_nanobody was called for each nanobody
    assert process_mock.call_count == 3
    process_mock.assert_any_call("nb1")
    process_mock.assert_any_call("nb2")
    process_mock.assert_any_call("nb3")

def test_process_skip_existing(dataset, monkeypatch):
    """Test process method skips already processed files"""
    # Mock _process_nanobody and torch.save
    process_mock = MagicMock(return_value=Data())
    monkeypatch.setattr(dataset, '_process_nanobody', process_mock)
    monkeypatch.setattr(torch, 'save', lambda data, path: None)

    # Mock os.path.exists to indicate files already exist
    monkeypatch.setattr(os.path, 'exists', lambda path: True)

    # Call process
    dataset.process()

    # Check that _process_nanobody was not called (files already exist)
    process_mock.assert_not_called()

def test_process_no_cache(dataset, monkeypatch):
    """Test process method is skipped when use_cache is False"""
    # Set use_cache to False
    dataset.use_cache = False

    # Mock _process_nanobody and torch.save
    process_mock = MagicMock(return_value=Data())
    monkeypatch.setattr(dataset, '_process_nanobody', process_mock)

    # Call process
    dataset.process()

    # Check that processing was skipped
    process_mock.assert_not_called()

def test_get_statistics(dataset, monkeypatch):
    """Test get_statistics method"""
    # Mock random choice for sampling
    monkeypatch.setattr(np.random, 'choice', lambda size, sample_size, replace: range(min(size, sample_size)))

    # Mock get method to return Data objects with different sizes
    data_objects = [
        Data(x=torch.zeros(10, 5), edge_index=torch.zeros(2, 15)),  # 10 nodes, 15 edges
        Data(x=torch.zeros(20, 5), edge_index=torch.zeros(2, 30)),  # 20 nodes, 30 edges
        Data(x=torch.zeros(15, 5), edge_index=torch.zeros(2, 25)),  # 15 nodes, 25 edges
    ]

    # Configure get method to return these data objects
    monkeypatch.setattr(dataset, 'get', lambda idx: data_objects[idx])

    # Get statistics
    stats = dataset.get_statistics()

    # Check statistics values
    assert stats["total_graphs"] == 3
    assert stats["avg_nodes"] == 15  # (10 + 20 + 15) / 3
    assert stats["min_nodes"] == 10
    assert stats["max_nodes"] == 20
    assert stats["avg_edges"] == 70 / 3  # (15 + 30 + 25) / 3
    assert stats["min_edges"] == 15
    assert stats["max_edges"] == 30

# Add this block only if running directly
if __name__ == "__main__":
    pytest.main(["-xvs", __file__])