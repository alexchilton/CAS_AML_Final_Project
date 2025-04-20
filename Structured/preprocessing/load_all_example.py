from ProteinGraphDataset import ProteinGraphDataset
import torch
from torch_geometric.loader import DataLoader
import matplotlib.pyplot as plt
import numpy as np
import time


def load_protein_dataset(root_path="/Users/alexchilton/Downloads/nanos_networkx_small"):
    """
    Load and prepare the protein graph dataset for GNN training.

    Args:
        root_path: Path to the directory containing protein data

    Returns:
        ProteinGraphDataset object
    """
    # Initialize dataset
    dataset = ProteinGraphDataset(
        root=root_path,
        use_cache=True
    )

    # Process all graphs (this will cache them for future use)
    print(f"Processing {len(dataset)} protein graphs...")
    start_time = time.time()
    dataset.process()
    processing_time = time.time() - start_time
    print(f"Processing completed in {processing_time:.2f} seconds")

    # Print dataset statistics
    stats = dataset.get_statistics()
    print("\nDataset Statistics:")
    for key, value in stats.items():
        if isinstance(value, (int, np.integer)):
            print(f"  {key}: {value}")
        else:
            print(f"  {key}: {value:.2f}")

    # Sample a few graphs
    print("\nSample graphs:")
    for i in range(min(3, len(dataset))):
        data = dataset[i]
        print(f"  Graph {i} - Nanobody: {data.nanobody_id}")
        print(f"    Nodes: {data.x.shape[0]}, Features: {data.x.shape[1]}")
        print(f"    Edges: {data.edge_index.shape[1]}, Edge Features: {data.edge_attr.shape[1]}")

    return dataset


def prepare_data_loaders(dataset, batch_size=32, split_ratio=[0.8, 0.1, 0.1]):
    """
    Prepare data loaders for training, validation, and testing.

    Args:
        dataset: ProteinGraphDataset object
        batch_size: Batch size for data loaders
        split_ratio: List of ratios for train/val/test split

    Returns:
        Dictionary containing train, val, and test data loaders
    """
    # Validate split ratio
    if sum(split_ratio) != 1.0:
        raise ValueError("Split ratios must sum to 1.0")

    # Calculate split sizes
    dataset_size = len(dataset)
    train_size = int(dataset_size * split_ratio[0])
    val_size = int(dataset_size * split_ratio[1])
    test_size = dataset_size - train_size - val_size

    # Create random indices for splitting
    indices = torch.randperm(dataset_size)
    train_indices = indices[:train_size]
    val_indices = indices[train_size:train_size+val_size]
    test_indices = indices[train_size+val_size:]

    # Create subset splits
    train_dataset = torch.utils.data.Subset(dataset, train_indices)
    val_dataset = torch.utils.data.Subset(dataset, val_indices)
    test_dataset = torch.utils.data.Subset(dataset, test_indices)

    # Create data loaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)

    # Print split information
    print(f"\nData split:")
    print(f"  Training: {len(train_dataset)} graphs")
    print(f"  Validation: {len(val_dataset)} graphs")
    print(f"  Testing: {len(test_dataset)} graphs")

    return {
        "train": train_loader,
        "val": val_loader,
        "test": test_loader
    }


def main():
    """Main function to demonstrate dataset loading and usage"""
    # Specify your data directory
    data_dir = "/Users/alexchilton/Downloads/nanos_networkx_small"

    # Load the dataset
    print(f"Loading protein graphs from: {data_dir}")
    dataset = load_protein_dataset(data_dir)

    # Validate and filter in one step
    filtered_dataset = validate_and_filter_dataset(dataset)


    # Validate the dataset
    print("\nValidating dataset consistency...")
    validation_results = validate_dataset_features(filtered_dataset)
    print_validation_results(validation_results)

    # Prepare data loaders for training
    data_loaders = prepare_data_loaders(filtered_dataset, batch_size=16)

    # Print a sample batch
    for batch in data_loaders["train"]:
        print("\nSample batch:")
        print(f"  Batch size: {batch.num_graphs}")
        print(f"  Node features: {batch.x.shape}")
        print(f"  Edge index: {batch.edge_index.shape}")
        print(f"  Edge attributes: {batch.edge_attr.shape}")
        # Only show one batch
        break
def validate_dataset_features(dataset, expected_node_features=39, expected_edge_features=2):
    """
    Validate that all graphs in the dataset have consistent and expected feature dimensions.

    Args:
        dataset: ProteinGraphDataset object
        expected_node_features: Expected number of node features (default: 39)
        expected_edge_features: Expected number of edge features (default: 2)

    Returns:
        Dictionary with validation results
    """
    problems = {
        "inconsistent_node_features": [],
        "inconsistent_edge_features": [],
        "missing_edge_attributes": [],
        "empty_graphs": []
    }

    for i in range(len(dataset)):
        data = dataset[i]
        nanobody_id = data.nanobody_id

        # Check for empty graphs
        if data.x.shape[0] == 0 or data.edge_index.shape[1] == 0:
            problems["empty_graphs"].append((i, nanobody_id))
            continue

        # Check node feature dimensions
        if data.x.shape[1] != expected_node_features:
            problems["inconsistent_node_features"].append((i, nanobody_id, data.x.shape[1]))

        # Check if edge attributes exist and have correct dimensions
        if not hasattr(data, 'edge_attr') or data.edge_attr.numel() == 0:
            problems["missing_edge_attributes"].append((i, nanobody_id))
        elif data.edge_attr.shape[0] > 0 and data.edge_attr.shape[1] != expected_edge_features:
            problems["inconsistent_edge_features"].append((i, nanobody_id, data.edge_attr.shape[1]))

    return problems
def validate_and_filter_dataset(dataset, expected_node_features=39, expected_edge_features=2):
    """
    Validates all graphs in the dataset and filters out any invalid ones.

    Args:
        dataset: ProteinGraphDataset object
        expected_node_features: Expected number of node features
        expected_edge_features: Expected number of edge features

    Returns:
        Filtered dataset without invalid graphs
    """
    import torch

    valid_indices = []
    invalid_ids = []

    # Check each graph
    for i in range(len(dataset)):
        data = dataset[i]
        nanobody_id = data.nanobody_id

        # Graph is valid if:
        # 1. It's not empty
        # 2. It has the correct node feature dimensions
        # 3. It has edge attributes with correct dimensions
        is_valid = (
                data.x.shape[0] > 0 and  # Not empty
                data.edge_index.shape[1] > 0 and  # Has edges
                data.x.shape[1] == expected_node_features and  # Correct node features
                hasattr(data, 'edge_attr') and  # Has edge attributes
                data.edge_attr.numel() > 0 and  # Edge attributes not empty
                data.edge_attr.shape[1] == expected_edge_features  # Correct edge features
        )

        if is_valid:
            valid_indices.append(i)
        else:
            invalid_ids.append(nanobody_id)

    # Create filtered dataset
    filtered_dataset = torch.utils.data.Subset(dataset, valid_indices)

    # Print results
    if invalid_ids:
        print(f"Removed {len(invalid_ids)} invalid graphs:")
        for id in invalid_ids:
            print(f"  - {id}")
    else:
        print("All graphs are valid!")

    print(f"Dataset size: {len(filtered_dataset)} (was {len(dataset)})")

    return filtered_dataset

# Example usage:
# dataset = load_protein_dataset(data_dir)
# filtered_dataset = validate_and_filter_dataset(dataset)
# data_loaders = prepare_data_loaders(filtered_dataset, batch_size=16)
def print_validation_results(validation_results):
    """
    Print a summary of dataset validation results.

    Args:
        validation_results: Dictionary with validation results
    """
    print("\nDataset Validation Results:")

    total_problems = sum(len(problems) for problems in validation_results.values())
    if total_problems == 0:
        print("  All graphs passed validation checks!")
        return

    print(f"  Found {total_problems} problems across the dataset:")

    if validation_results["empty_graphs"]:
        print(f"  - Empty graphs: {len(validation_results['empty_graphs'])}")
        for idx, nanobody_id in validation_results["empty_graphs"][:5]:
            print(f"    * Graph {idx} (Nanobody: {nanobody_id})")
        if len(validation_results["empty_graphs"]) > 5:
            print(f"      ... and {len(validation_results['empty_graphs']) - 5} more")

    if validation_results["inconsistent_node_features"]:
        print(f"  - Inconsistent node features: {len(validation_results['inconsistent_node_features'])}")
        for idx, nanobody_id, dim in validation_results["inconsistent_node_features"][:5]:
            print(f"    * Graph {idx} (Nanobody: {nanobody_id}) has {dim} features")
        if len(validation_results["inconsistent_node_features"]) > 5:
            print(f"      ... and {len(validation_results['inconsistent_node_features']) - 5} more")

    if validation_results["missing_edge_attributes"]:
        print(f"  - Missing edge attributes: {len(validation_results['missing_edge_attributes'])}")
        for idx, nanobody_id in validation_results["missing_edge_attributes"][:5]:
            print(f"    * Graph {idx} (Nanobody: {nanobody_id})")
        if len(validation_results["missing_edge_attributes"]) > 5:
            print(f"      ... and {len(validation_results['missing_edge_attributes']) - 5} more")

    if validation_results["inconsistent_edge_features"]:
        print(f"  - Inconsistent edge features: {len(validation_results['inconsistent_edge_features'])}")
        for idx, nanobody_id, dim in validation_results["inconsistent_edge_features"][:5]:
            print(f"    * Graph {idx} (Nanobody: {nanobody_id}) has {dim} features")
        if len(validation_results["inconsistent_edge_features"]) > 5:
            print(f"      ... and {len(validation_results['inconsistent_edge_features']) - 5} more")

if __name__ == "__main__":
    main()