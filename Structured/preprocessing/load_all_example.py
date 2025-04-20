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

    # Prepare data loaders for training
    data_loaders = prepare_data_loaders(dataset, batch_size=16)

    # Print a sample batch
    for batch in data_loaders["train"]:
        print("\nSample batch:")
        print(f"  Batch size: {batch.num_graphs}")
        print(f"  Node features: {batch.x.shape}")
        print(f"  Edge index: {batch.edge_index.shape}")
        print(f"  Edge attributes: {batch.edge_attr.shape}")
        # Only show one batch
        break


if __name__ == "__main__":
    main()