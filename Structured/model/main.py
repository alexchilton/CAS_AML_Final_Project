import torch
import numpy as np
import os
import matplotlib.pyplot as plt
from DatasetHandler import DatasetHandler
from GATVAE import GATVAE
from TrainingManager import TrainingManager
from LossVisualizer import LossVisualizer
from LatentSpaceManager import LatentSpaceManager
from ModelUtils import ModelUtils
from torch_geometric.loader import DataLoader as PyGDataLoader


def main():
    # Configuration
    data_root = "/Users/alexchilton/Downloads/nanos_networkx_small"  # Path to your dataset
    batch_size = 16  # Smaller batch size to avoid memory issues
    hidden_dim = 64
    latent_dim = 32
    num_heads = 4
    learning_rate = 0.001
    num_epochs = 30
    kl_weight = 0.01  # Reduced KL weight to prevent vanishing KL loss

    # Set random seed for reproducibility
    torch.manual_seed(42)
    np.random.seed(42)

    # Create output directories
    os.makedirs("./models", exist_ok=True)
    os.makedirs("./plots", exist_ok=True)
    os.makedirs("./latent_plots", exist_ok=True)

    # Initialize dataset handler
    dataset_handler = DatasetHandler(data_root=data_root, batch_size=batch_size)

    # Load and prepare dataset
    # Load dataset with global normalization
    dataset = dataset_handler.load_dataset(use_global_normalization=True)

    # Save parameters for future use
    dataset_handler.save_global_normalization_params("params.json")
    train_loader, val_loader, test_loader = dataset_handler.prepare_dataloaders()

    # Get input dimension
    input_dim = dataset_handler.get_input_dim()
    print(f"Input dimension: {input_dim}")

    # Create model
    model = GATVAE(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        latent_dim=latent_dim,
        num_heads=num_heads
    )

    # Print model summary
    print(model)
    num_params = sum(p.numel() for p in model.parameters())
    print(f"Number of parameters: {num_params}")

    # Create training manager
    trainer = TrainingManager(
        model=model,
        learning_rate=learning_rate,
        kl_weight=kl_weight
    )

    # Create model save path
    model_save_path = "./models/gatvae_model.pth"

    # Train model
    print("Training model...")
    history = trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        num_epochs=num_epochs,
        save_path=model_save_path
    )

    # Create loss visualizer
    visualizer = LossVisualizer(save_dir="./plots")

    # Plot losses
    visualizer.plot_all(history)

    # Save model with configuration
    config = {
        "input_dim": input_dim,
        "hidden_dim": hidden_dim,
        "latent_dim": latent_dim,
        "num_heads": num_heads,
        "learning_rate": learning_rate,
        "kl_weight": kl_weight,
        "num_epochs": num_epochs
    }

    ModelUtils.save_model(model, model_save_path, config)

    # Analyze latent space
    latent_manager = LatentSpaceManager(model, device=trainer.device, save_dir="./latent_plots")

    # Use a smaller subset for visualization to avoid memory issues
    subset_size = min(100, len(test_loader.dataset))
    subset_indices = np.random.choice(len(test_loader.dataset), subset_size, replace=False)
    subset_dataset = torch.utils.data.Subset(test_loader.dataset, subset_indices)

    # Use PyTorch Geometric's DataLoader instead of PyTorch's standard DataLoader
    subset_loader = PyGDataLoader(subset_dataset, batch_size=8)

    # Encode test subset
    print("Encoding protein subset for latent space visualization...")
    latent_vectors, ids = latent_manager.encode_dataset(subset_loader)

    # Visualize latent space
    print("Visualizing latent space...")
    latent_manager.visualize_2d(latent_vectors, method='pca')
    latent_manager.visualize_2d(latent_vectors, method='tsne')

    # Export latent vectors
    ModelUtils.export_latent_vectors(latent_vectors, ids, "./latent_plots/latent_vectors.npz")

    # Sample from latent space (if dataset not empty)
    try:
        print("Sampling from latent space...")
        if len(dataset) > 0:
            sample_data = dataset[0].to(trainer.device)
            sampled_graphs = latent_manager.sample_from_latent(
                num_samples=5,
                template_data=sample_data
            )
            print(f"Generated {len(sampled_graphs)} samples from latent space")
    except Exception as e:
        print(f"Error sampling from latent space: {e}")

    # Try to interpolate between two proteins
    try:
        print("Interpolating between proteins...")
        if len(dataset) >= 2:
            protein1 = dataset[0].to(trainer.device)
            protein2 = dataset[1].to(trainer.device)
            interpolated = latent_manager.interpolate(protein1, protein2, steps=5)
            print(f"Generated {len(interpolated)} interpolated proteins")
    except Exception as e:
        print(f"Error interpolating proteins: {e}")

    print("Done!")


if __name__ == "__main__":

    main()