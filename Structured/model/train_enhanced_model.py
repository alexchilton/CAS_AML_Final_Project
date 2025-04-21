import torch
import numpy as np
import os
from torch_geometric.loader import DataLoader as PyGDataLoader

# Import dataset handling
from DatasetHandler import DatasetHandler

# Import enhanced model components
from EnhancedGATEncoder import EnhancedGATEncoder
from EnhancedGATDecoder import EnhancedGATDecoder
from EnhancedGATVAE import EnhancedGATVAE
from EnhancedTrainingManager import EnhancedTrainingManager


import os
import argparse
import matplotlib.pyplot as plt
import sys  # You have sys imported here

# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../generation'))
# Import visualization components
from LossVisualizer import LossVisualizer
from LatentSpaceManager import LatentSpaceManager
from GenerationVisualizer import GenerationVisualizer


def main():
    """Main function to demonstrate training the enhanced model."""
    # Configuration
    data_root = "/Users/alexchilton/Downloads/nanos_networkx_small"  # Update this path
    batch_size = 16
    hidden_dim = 64
    latent_dim = 32
    num_heads = 4
    pos_embedding_dim = 16
    spatial_embedding_dim = 16
    edge_hidden_dim = 64
    learning_rate = 0.001
    num_epochs = 30
    kl_weight = 0.01
    edge_weight = 0.5
    spatial_weight = 0.2
    lambda_reg = 0.001

    # Set random seed for reproducibility
    torch.manual_seed(42)
    np.random.seed(42)

    # Create output directories
    os.makedirs("./models", exist_ok=True)
    os.makedirs("./plots", exist_ok=True)
    os.makedirs("./latent_plots", exist_ok=True)
    os.makedirs("./generation_visualizations", exist_ok=True)

    # Initialize dataset handler
    print(f"Loading dataset from {data_root}")
    dataset_handler = DatasetHandler(data_root=data_root, batch_size=batch_size)

    # Load and prepare dataset
    dataset = dataset_handler.load_dataset()
    train_loader, val_loader, test_loader = dataset_handler.prepare_dataloaders()

    # Get input dimension
    input_dim = dataset_handler.get_input_dim()
    print(f"Input dimension: {input_dim}")

    # Create enhanced model
    model = EnhancedGATVAE(
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        latent_dim=latent_dim,
        num_heads=num_heads,
        pos_embedding_dim=pos_embedding_dim,
        spatial_embedding_dim=spatial_embedding_dim,
        edge_hidden_dim=edge_hidden_dim
    )

    # Print model summary
    print(model)
    num_params = sum(p.numel() for p in model.parameters())
    print(f"Number of parameters: {num_params}")

    # Create enhanced training manager
    trainer = EnhancedTrainingManager(
        model=model,
        learning_rate=learning_rate,
        kl_weight=kl_weight,
        edge_weight=edge_weight,
        spatial_weight=spatial_weight,
        lambda_reg=lambda_reg
    )

    # Create model save path
    model_save_path = "./models/enhanced_gatvae_model.pth"

    # Train model
    print("Training model...")
    history = trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        num_epochs=num_epochs,
        save_path=model_save_path
    )

    # Create loss visualizer and plot training history
    print("Visualizing training history...")
    loss_visualizer = LossVisualizer(save_dir="./plots")

    # Plot training and validation loss
    loss_visualizer.plot_training_history({
        'train_losses': history['train_loss'],
        'val_losses': history['val_loss']
    })

    # Plot component losses
    loss_visualizer.plot_component_losses({
        'reconstruction_losses': history['reconstruction_loss'],
        'kl_losses': history['kl_loss'],
        'edge_losses': history['edge_loss'],
        'spatial_losses': history['spatial_loss']
    })

    # Plot loss ratios
    loss_visualizer.plot_loss_ratio({
        'reconstruction_losses': history['reconstruction_loss'],
        'kl_losses': history['kl_loss']
    })

    # Evaluate on test data
    print("Evaluating model on test data...")
    metrics = trainer.evaluate(test_loader)
    print("Test metrics:")
    for key, value in metrics.items():
        print(f"  {key}: {value:.4f}")

    # Initialize latent space manager
    latent_manager = LatentSpaceManager(
        model=model,
        device=trainer.device,
        save_dir="./latent_plots"
    )

    # Visualize latent space
    print("Visualizing latent space...")
    # Use a smaller subset for visualization
    subset_size = min(100, len(test_loader.dataset))
    subset_indices = np.random.choice(len(test_loader.dataset), subset_size, replace=False)
    subset_dataset = torch.utils.data.Subset(test_loader.dataset, subset_indices)
    subset_loader = PyGDataLoader(subset_dataset, batch_size=8)

    # Encode test subset
    latent_vectors, ids = latent_manager.encode_dataset(subset_loader)

    # Visualize latent space with PCA and t-SNE
    latent_manager.visualize_2d(latent_vectors, method='pca')
    latent_manager.visualize_2d(latent_vectors, method='tsne')

    # Initialize generation visualizer
    generation_visualizer = GenerationVisualizer(
        model=model,
        device=trainer.device,
        save_dir="./generation_visualizations"
    )

    # Generate new proteins
    print("Generating new protein structures...")
    generated_proteins = model.sample(num_samples=5, num_nodes=20, device=trainer.device)

    # Visualize generated proteins
    for i, protein in enumerate(generated_proteins):
        generation_visualizer.plot_protein_graph(
            protein,
            title=f"Generated Protein {i+1}",
            save=True,
            filename=f"generated_protein_{i+1}.png"
        )

    # Interpolate between proteins if we have enough test proteins
    if len(test_loader.dataset) >= 2:
        print("Generating protein interpolation...")

        # Get two proteins from test set
        protein1 = test_loader.dataset[0].to(trainer.device)
        protein2 = test_loader.dataset[1].to(trainer.device)

        # Generate interpolation
        interpolated_proteins = latent_manager.interpolate(protein1, protein2, steps=8)

        # Visualize interpolation
        generation_visualizer.plot_interpolation_sequence(
            interpolated_proteins,
            save=True,
            filename="protein_interpolation.png"
        )

        print(f"Generated {len(interpolated_proteins)} interpolated proteins")

    print("Done!")


if __name__ == "__main__":
    main()