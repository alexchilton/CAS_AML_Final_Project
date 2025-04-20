import torch
import torch.nn.functional as F
import numpy as np
from typing import List, Dict, Any, Optional, Tuple, Callable
from torch_geometric.data import Data
from tqdm import tqdm
import os
import sys

# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../preprocessing'))

from GATVAE import GATVAE


class GenerationOptimizer:
    """
    Class for optimizing protein structures in the latent space.
    Provides methods for latent space exploration and optimization.
    """

    def __init__(self, model: GATVAE, device='cpu'):
        """
        Initialize the generation optimizer.

        Args:
            model: Trained GATVAE model
            device: Device to run computations on
        """
        self.model = model
        self.device = device

        # Move model to device and set to evaluation mode
        self.model.to(device)
        self.model.eval()

    def optimize_latent_vector(self, initial_z: torch.Tensor, objective_fn: Callable,
                               template_data: Data, lr: float = 0.01, steps: int = 100,
                               verbose: bool = True) -> Tuple[torch.Tensor, List[float]]:
        """
        Optimize a latent vector based on an objective function.

        Args:
            initial_z: Initial latent vector to optimize
            objective_fn: Function that evaluates protein quality (higher is better)
            template_data: Template graph structure to use for reconstruction
            lr: Learning rate for optimization
            steps: Number of optimization steps
            verbose: Whether to show progress

        Returns:
            Tuple of (optimized_latent_vector, objective_values)
        """
        # Create a copy of initial_z that requires gradients
        z = initial_z.clone().to(self.device).requires_grad_(True)
        template_data = template_data.to(self.device)

        # Optimizer
        optimizer = torch.optim.Adam([z], lr=lr)

        # Track objective values
        objective_values = []

        # Optimization loop
        iterator = tqdm(range(steps)) if verbose else range(steps)
        for i in iterator:
            # Reset gradients
            optimizer.zero_grad()

            # Decode current latent vector
            with torch.set_grad_enabled(True):
                x_reconstructed = self.model.decode(z.unsqueeze(0), template_data)

                # Create protein graph
                protein = Data(
                    x=x_reconstructed,
                    edge_index=template_data.edge_index,
                    edge_attr=template_data.edge_attr
                )

                # Calculate objective (negative because we're minimizing)
                obj_value = objective_fn(protein)
                objective_values.append(obj_value.item())

                # We want to maximize the objective, so minimize negative objective
                loss = -obj_value

                # Backpropagate
                loss.backward()

                # Update latent vector
                optimizer.step()

            if verbose and (i+1) % 10 == 0:
                tqdm.write(f"Step {i+1}, Objective: {obj_value.item():.4f}")

        # Return optimized latent vector (detached) and objective history
        return z.detach(), objective_values

    def latent_space_exploration(self, base_z: torch.Tensor, template_data: Data,
                                 num_directions: int = 10, step_size: float = 0.5,
                                 steps_per_direction: int = 5) -> List[Data]:
        """
        Explore the latent space by moving in random directions from a base point.

        Args:
            base_z: Base latent vector to start from
            template_data: Template graph structure to use for reconstruction
            num_directions: Number of random directions to explore
            step_size: Size of steps in each direction
            steps_per_direction: Number of steps to take in each direction

        Returns:
            List of protein graphs from exploration
        """
        self.model.eval()
        base_z = base_z.to(self.device)
        template_data = template_data.to(self.device)

        # Generate random unit directions
        directions = []
        for _ in range(num_directions):
            direction = torch.randn_like(base_z)
            direction = direction / direction.norm()  # Normalize to unit vector
            directions.append(direction)

        # Explore each direction
        explored_proteins = []

        with torch.no_grad():
            # First add the base protein
            x_reconstructed = self.model.decode(base_z.unsqueeze(0), template_data)
            base_protein = Data(
                x=x_reconstructed,
                edge_index=template_data.edge_index,
                edge_attr=template_data.edge_attr
            )
            explored_proteins.append(base_protein)

            # Explore each direction
            for direction in directions:
                for step in range(1, steps_per_direction + 1):
                    # Move in the direction
                    current_z = base_z + direction * (step * step_size)

                    # Decode
                    x_reconstructed = self.model.decode(current_z.unsqueeze(0), template_data)

                    # Create protein
                    protein = Data(
                        x=x_reconstructed,
                        edge_index=template_data.edge_index,
                        edge_attr=template_data.edge_attr
                    )
                    explored_proteins.append(protein)

        return explored_proteins

    def feature_guided_generation(self, target_features: torch.Tensor, template_data: Data,
                                  lr: float = 0.01, steps: int = 100,
                                  weight_decay: float = 0.001) -> Data:
        """
        Generate a protein with features similar to the target features.

        Args:
            target_features: Target features to guide generation
            template_data: Template graph structure to use for reconstruction
            lr: Learning rate for optimization
            steps: Number of optimization steps
            weight_decay: Regularization parameter

        Returns:
            Generated protein graph
        """
        self.model.eval()
        template_data = template_data.to(self.device)
        target_features = target_features.to(self.device)

        # Initialize random latent vector
        z = torch.randn(self.model.latent_dim, device=self.device, requires_grad=True)

        # Optimizer with weight decay for regularization
        optimizer = torch.optim.Adam([z], lr=lr, weight_decay=weight_decay)

        # Optimization loop
        for step in tqdm(range(steps)):
            # Reset gradients
            optimizer.zero_grad()

            # Decode current latent vector
            with torch.set_grad_enabled(True):
                x_reconstructed = self.model.decode(z.unsqueeze(0), template_data)

                # Calculate feature loss (MSE between generated and target features)
                # Note: We may need to average/aggregate node features depending on the target
                feature_loss = F.mse_loss(x_reconstructed, target_features)

                # Backpropagate
                feature_loss.backward()

                # Update latent vector
                optimizer.step()

            if (step+1) % 10 == 0:
                print(f"Step {step+1}, Loss: {feature_loss.item():.4f}")

        # Generate final protein using optimized latent vector
        with torch.no_grad():
            x_reconstructed = self.model.decode(z.unsqueeze(0), template_data)
            generated_protein = Data(
                x=x_reconstructed,
                edge_index=template_data.edge_index,
                edge_attr=template_data.edge_attr
            )

        return generated_protein

    def ensemble_latent_vectors(self, latent_vectors: List[torch.Tensor],
                                weights: Optional[List[float]] = None) -> torch.Tensor:
        """
        Create an ensemble of latent vectors with optional weighting.

        Args:
            latent_vectors: List of latent vectors to ensemble
            weights: Optional list of weights for each vector (default: equal weights)

        Returns:
            Ensembled latent vector
        """
        # Convert list to tensor
        stacked_vectors = torch.stack([z.to(self.device) for z in latent_vectors])

        # Apply weights if provided, otherwise use equal weights
        if weights is None:
            ensembled_z = torch.mean(stacked_vectors, dim=0)
        else:
            # Normalize weights to sum to 1
            weights = torch.tensor(weights, device=self.device)
            weights = weights / weights.sum()

            # Apply weights
            ensembled_z = torch.sum(stacked_vectors * weights.view(-1, 1), dim=0)

        return ensembled_z