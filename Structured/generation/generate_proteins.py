import torch
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
import sys  # You have sys imported here

# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../model'))

from DatasetHandler import DatasetHandler
from GATVAE import GATVAE
from ModelUtils import ModelUtils
from torch_geometric.loader import DataLoader as PyGDataLoader

# Import the new generation classes
from ProteinGenerator import ProteinGenerator
from GenerationOptimizer import GenerationOptimizer
from GenerationVisualizer import GenerationVisualizer


def parse_args():
    parser = argparse.ArgumentParser(description='Protein Graph Generation')
    parser.add_argument('--data_dir', type=str, required=True,
                        help='Path to the dataset directory')
    parser.add_argument('--model_path', type=str, default='../model/models/gatvae_model.pth',
                        help='Path to the trained model file')
    parser.add_argument('--output_dir', type=str, default='./generation_results',
                        help='Output directory for generated results')
    parser.add_argument('--num_samples', type=int, default=5,
                        help='Number of proteins to generate')
    parser.add_argument('--mode', type=str, default='random',
                        choices=['random', 'interpolate', 'optimize', 'explore', 'average'],
                        help='Generation mode')
    parser.add_argument('--batch_size', type=int, default=16,
                        help='Batch size for data loading')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    return parser.parse_args()


def main():
    # Parse command line arguments
    args = parse_args()

    # Set random seed for reproducibility
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Set device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Load dataset
    print(f"Loading protein dataset from {args.data_dir}")
    dataset_handler = DatasetHandler(data_root=args.data_dir, batch_size=args.batch_size)
    dataset = dataset_handler.load_dataset()
    train_loader, val_loader, test_loader = dataset_handler.prepare_dataloaders()


    input_dim = dataset_handler.get_input_dim()

    # Load model
    print(f"Loading model from {args.model_path}")
    model = ModelUtils.load_model(args.model_path, input_dim)
    model.to(device)
    model.eval()

    # Initialize generation classes
    protein_generator = ProteinGenerator(model, device=device)
    generation_optimizer = GenerationOptimizer(model, device=device)
    generation_visualizer = GenerationVisualizer(
        model,
        device=device,
        save_dir=args.output_dir
    )

    # Select proteins to use as templates/sources
    test_loader = dataset_handler.test_loader
    test_proteins = []
    for i in range(min(args.num_samples, len(test_loader.dataset))):
        test_proteins.append(test_loader.dataset[i])

    if not test_proteins:
        print("Error: No proteins available in the dataset")
        return

    # Generate proteins based on the selected mode
    if args.mode == 'random':
        generate_random_proteins(
            protein_generator,
            generation_visualizer,
            test_proteins,
            args.num_samples,
            device
        )
    elif args.mode == 'interpolate':
        generate_interpolated_proteins(
            protein_generator,
            generation_visualizer,
            test_proteins,
            device
        )
    elif args.mode == 'optimize':
        generate_optimized_proteins(
            protein_generator,
            generation_optimizer,
            generation_visualizer,
            model,
            test_proteins,
            device,
            latent_dim=model.latent_dim
        )
    elif args.mode == 'explore':
        explore_latent_space(
            protein_generator,
            generation_optimizer,
            generation_visualizer,
            model,
            test_proteins,
            device
        )
    elif args.mode == 'average':
        average_proteins(
            protein_generator,
            generation_visualizer,
            test_proteins,
            device
        )

    print(f"Protein generation completed. Results saved to {args.output_dir}")


def generate_random_proteins(generator, visualizer, test_proteins, num_samples, device):
    print(f"Generating {num_samples} random protein samples...")
    template_protein = test_proteins[0].to(device)

    # Generate random protein samples
    random_samples = generator.sample_random(
        num_samples=num_samples,
        template_data=template_protein
    )

    # Visualize the generated proteins
    for i, protein in enumerate(random_samples):
        visualizer.plot_protein_graph(
            protein,
            title=f"Random Protein {i+1}",
            save=True,
            filename=f"random_protein_{i+1}.png"
        )

    print(f"Generated {len(random_samples)} random protein samples")


def generate_interpolated_proteins(generator, visualizer, test_proteins, device):
    if len(test_proteins) < 2:
        print("Error: Need at least 2 proteins for interpolation")
        return

    print("Generating interpolated proteins...")
    protein1 = test_proteins[0].to(device)
    protein2 = test_proteins[1].to(device)

    # Generate interpolated proteins
    steps = 8
    interpolated = generator.interpolate(protein1, protein2, steps=steps)

    # Visualize individual proteins
    for i, protein in enumerate(interpolated):
        visualizer.plot_protein_graph(
            protein,
            title=f"Interpolation Step {i+1}/{steps}",
            show=False,
            save=True,
            filename=f"interpolation_step_{i+1}.png"
        )

    # Visualize the entire sequence
    visualizer.plot_interpolation_sequence(
        interpolated,
        save=True,
        filename="interpolation_sequence.png"
    )

    print(f"Generated {len(interpolated)} interpolated proteins")


def generate_optimized_proteins(generator, optimizer, visualizer, model, test_proteins, device, latent_dim):
    print("Generating optimized proteins...")
    template_protein = test_proteins[0].to(device)

    # Define a simple objective function
    def objective_fn(protein):
        # Example: reward higher average feature values in certain dimensions
        # In a real application, this would be a meaningful protein quality metric
        target_feature_idx = 25  # Example feature index
        feature_mean = protein.x[:, target_feature_idx].mean()
        return feature_mean

    # Create initial random latent vector
    initial_z = torch.randn(latent_dim, device=device)

    # Optimize the latent vector
    optimized_z, obj_values = optimizer.optimize_latent_vector(
        initial_z=initial_z,
        objective_fn=objective_fn,
        template_data=template_protein,
        lr=0.01,
        steps=100,
        verbose=True
    )

    # Generate proteins using initial and optimized vectors
    with torch.no_grad():
        initial_protein = generator.reconstruct(template_protein)

        optimized_protein_z = model.decode(optimized_z.unsqueeze(0), template_protein)
        optimized_protein = torch.utils.data.dataset.Data(
            x=optimized_protein_z,
            edge_index=template_protein.edge_index,
            edge_attr=template_protein.edge_attr
        )

    # Visualize the original, initial, and optimized proteins
    visualizer.compare_original_vs_generated(
        template_protein,
        optimized_protein,
        save=True,
        filename="original_vs_optimized.png"
    )

    # Plot optimization progress
    plt.figure(figsize=(10, 6))
    plt.plot(obj_values)
    plt.title('Optimization Progress')
    plt.xlabel('Step')
    plt.ylabel('Objective Value')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(visualizer.save_dir, 'optimization_curve.png'))
    plt.close()

    print("Optimization completed")


def explore_latent_space(generator, optimizer, visualizer, model, test_proteins, device):
    print("Exploring latent space neighborhood...")
    template_protein = test_proteins[0].to(device)

    # Encode a protein to get a starting point
    with torch.no_grad():
        base_z = model.encode(template_protein)

    # Explore the latent space neighborhood
    explored_proteins = optimizer.latent_space_exploration(
        base_z=base_z,
        template_data=template_protein,
        num_directions=5,
        step_size=0.5,
        steps_per_direction=4
    )

    # Visualize some of the explored proteins
    for i, protein in enumerate(explored_proteins[:10]):
        visualizer.plot_protein_graph(
            protein,
            title=f"Exploration Sample {i+1}",
            show=False,
            save=True,
            filename=f"exploration_sample_{i+1}.png"
        )

    print(f"Generated {len(explored_proteins)} proteins through latent space exploration")


def average_proteins(generator, visualizer, test_proteins, device):
    if len(test_proteins) < 2:
        print("Error: Need at least 2 proteins for averaging")
        return

    num_to_average = min(3, len(test_proteins))
    print(f"Averaging {num_to_average} proteins in latent space...")

    proteins_to_average = test_proteins[:num_to_average]

    # Generate the averaged protein
    averaged_protein = generator.average_proteins(proteins_to_average)

    # Visualize individual proteins that were averaged
    for i, protein in enumerate(proteins_to_average):
        visualizer.plot_protein_graph(
            protein,
            title=f"Source Protein {i+1}",
            show=False,
            save=True,
            filename=f"source_protein_{i+1}.png"
        )

    # Visualize the averaged protein
    visualizer.plot_protein_graph(
        averaged_protein,
        title="Averaged Protein",
        save=True,
        filename="averaged_protein.png"
    )

    print("Protein averaging completed")


if __name__ == "__main__":
    main()