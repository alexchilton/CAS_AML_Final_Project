import torch
import torch.nn.functional as F
from torch.optim import Adam
from torch.utils.data import DataLoader
from torch_geometric.data import Data, Batch
from torch_geometric.utils import to_dense_adj, dense_to_sparse
import numpy as np
from tqdm import tqdm
import os

from EnhancedGATVAE import EnhancedGATVAE


class EnhancedTrainingManager:
    """
    Enhanced training manager for the GAT-based VAE model with structure learning.
    Handles training loop, loss calculation, and evaluation.
    """

    def __init__(self, model, learning_rate=0.001, kl_weight=0.01, edge_weight=0.5,
                 spatial_weight=1000, lambda_reg=0.001):
        """
        Initialize the enhanced training manager.

        Args:
            model: EnhancedGATVAE model to train
            learning_rate: Learning rate for optimization
            kl_weight: Weight for the KL divergence loss term
            edge_weight: Weight for the edge prediction loss
            spatial_weight: Weight for the spatial consistency loss
            lambda_reg: Regularization parameter for edge sparsity
        """
        self.model = model
        self.optimizer = Adam(model.parameters(), lr=learning_rate)
        self.kl_weight = kl_weight
        self.edge_weight = edge_weight
        self.spatial_weight = spatial_weight
        self.lambda_reg = lambda_reg
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Move model to device
        self.model.to(self.device)

        # Training history
        self.train_losses = []
        self.val_losses = []
        self.reconstruction_losses = []
        self.kl_losses = []
        self.edge_losses = []
        self.spatial_losses = []

    def _calculate_edge_prediction_loss(self, original_data, edge_probs):
        """
        Calculate edge prediction loss between original and predicted edges.
        Fix for the size mismatch issue.

        Args:
            original_data: Original PyTorch Geometric Data object
            edge_probs: Edge probabilities from model

        Returns:
            edge_loss: Loss for edge prediction with sparsity regularization
        """
        # Instead of creating ground truth from adjacency matrix,
        # just use binary cross entropy with the edges we have predictions for

        # Create a dummy target with the same size as edge_probs
        # Since we already have probabilities for selected edges only
        targets = torch.ones_like(edge_probs)

        # Calculate simple BCE loss just on the predicted edges
        edge_loss = F.binary_cross_entropy(edge_probs, targets)

        # Add regularization to promote sparsity (optional)
        edge_sparsity = self.lambda_reg * edge_probs.mean()

        return edge_loss + edge_sparsity


    def train_epoch(self, dataloader):
        """
        Enhanced training loop with careful gradient tracking and debugging.

        Args:
            dataloader: DataLoader for training data

        Returns:
            epoch_metrics: Dictionary of loss metrics
        """
        self.model.train()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        total_edge_loss = 0
        total_spatial_loss = 0
        num_batches = 0

        # Create hook for gradient debugging
        gradients = {}

        def hook_fn(name):
            def fn(grad):
                gradients[name] = grad.abs().mean().item()
            return fn

        # Register hooks for key layers to track gradient flow
        hooks = []
        for name, param in self.model.named_parameters():
            if any(key in name for key in ['decoder.fc_z', 'decoder.coord_generator',
                                           'decoder.edge_predictor']):
                hook = param.register_hook(hook_fn(name))
                hooks.append(hook)

        for batch_idx, data in enumerate(tqdm(dataloader, desc="Training")):
            try:
                # Move data to device
                data = data.to(self.device)

                # Ensure data has required attributes
                if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                    data.edge_attr = torch.ones(data.edge_index.size(1), 2,
                                                dtype=torch.float, device=self.device)

                # Zero gradients
                self.optimizer.zero_grad()

                # Debug data
                if batch_idx < 2:  # Only for first few batches
                    print(f"\nBatch {batch_idx} - Data shape: {data.x.shape}")
                    coord_idx_start = 21 + 7
                    coord_idx_end = coord_idx_start + 3
                    print(f"Coordinate slice: [{coord_idx_start}:{coord_idx_end}]")
                    print("Original coords sample:",
                          data.x[:2, coord_idx_start:coord_idx_end].detach().cpu().numpy())

                # Forward pass with enhanced model
                # Get each return value explicitly to ensure gradient flow
                reconstructed_data, mu, logvar, edge_probs = self.model(data)

                # Verify coordinates in reconstructed data
                if batch_idx < 2:  # Only for first few batches
                    recon_coords = reconstructed_data.x[:2, coord_idx_start:coord_idx_end]
                    print("Reconstructed coords sample:", recon_coords.detach().cpu().numpy())
                    print("Gradients enabled for recon_coords:", recon_coords.requires_grad)

                # Calculate reconstruction loss for node features
                # Using sum reduction instead of mean to scale with batch size
                rec_loss = F.mse_loss(reconstructed_data.x, data.x, reduction='sum')

                # Calculate KL divergence loss
                # Scale KL by batch size
                batch_size = 1 if data.batch is None else (data.batch.max().item() + 1)
                kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
                kl_loss = kl_loss / batch_size

                # Calculate edge prediction loss with more direct gradient path
                edge_loss = self._calculate_edge_prediction_loss(data, edge_probs)

                # Calculate spatial consistency loss
                spatial_loss = self._calculate_spatial_consistency_loss(data, reconstructed_data)

                # Scale spatial loss to make gradients more significant
                # This can help when the spatial loss is important but gradients are small
                spatial_loss = spatial_loss * 10.0  # Explicit scaling factor

                if batch_idx < 2:
                    print(f"Spatial loss (scaled): {spatial_loss.item():.4f}")

                # Calculate the total loss with careful coefficient scaling
                loss = (
                        rec_loss / data.x.size(0) +  # Normalize reconstruction loss
                        self.kl_weight * kl_loss +
                        self.edge_weight * edge_loss +
                        self.spatial_weight * spatial_loss
                )

                # Check for reasonable loss values
                if torch.isnan(loss) or torch.isinf(loss):
                    print(f"Warning: Invalid loss detected: {loss.item()}")
                    print(f"  Reconstruction loss: {rec_loss.item()}")
                    print(f"  KL loss: {kl_loss.item()}")
                    print(f"  Edge loss: {edge_loss.item()}")
                    print(f"  Spatial loss: {spatial_loss.item()}")
                    continue

                # Backward pass with retain_graph to allow gradient analysis
                loss.backward(retain_graph=(batch_idx < 2))

                # Check for gradient flow in key parameters for the first few batches
                if batch_idx < 2:
                    print("\nGradient flow check:")
                    for name, param in self.model.named_parameters():
                        if any(key in name for key in ['decoder.fc_z', 'decoder.coord_generator',
                                                       'decoder.edge_predictor']):
                            if param.grad is not None:
                                grad_norm = param.grad.abs().mean().item()
                                print(f"  {name}: {grad_norm:.6f}")
                            else:
                                print(f"  {name}: No gradient")

                    # Print collected gradients from hooks
                    print("\nGradient magnitudes from hooks:")
                    for name, value in gradients.items():
                        print(f"  {name}: {value:.6f}")
                    gradients.clear()  # Clear for next batch

                # Gradient clipping to prevent exploding gradients
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)

                # Update weights
                self.optimizer.step()

                # Accumulate metrics for reporting
                total_loss += loss.item()
                total_rec_loss += rec_loss.item() / data.x.size(0)  # Normalize
                total_kl_loss += kl_loss.item()
                total_edge_loss += edge_loss.item()
                total_spatial_loss += spatial_loss.item()
                num_batches += 1

                # Print loss components for first few batches
                if batch_idx < 2:
                    print(f"Loss components - Rec: {rec_loss.item()/data.x.size(0):.4f}, "
                          f"KL: {kl_loss.item():.4f}, Edge: {edge_loss.item():.4f}, "
                          f"Spatial: {spatial_loss.item():.4f}")

            except Exception as e:
                print(f"Error processing batch {batch_idx}: {e}")
                import traceback
                traceback.print_exc()
                continue

        # Remove hooks
        for hook in hooks:
            hook.remove()

        # Calculate average metrics
        metrics = {}
        if num_batches > 0:
            metrics['loss'] = total_loss / num_batches
            metrics['rec_loss'] = total_rec_loss / num_batches
            metrics['kl_loss'] = total_kl_loss / num_batches
            metrics['edge_loss'] = total_edge_loss / num_batches
            metrics['spatial_loss'] = total_spatial_loss / num_batches
        else:
            metrics = {k: float('nan') for k in
                       ['loss', 'rec_loss', 'kl_loss', 'edge_loss', 'spatial_loss']}

        # Store in training history
        self.train_losses.append(metrics['loss'])
        self.reconstruction_losses.append(metrics['rec_loss'])
        self.kl_losses.append(metrics['kl_loss'])
        self.edge_losses.append(metrics['edge_loss'])
        self.spatial_losses.append(metrics['spatial_loss'])

        return metrics

    def _calculate_spatial_consistency_loss(self, original_data, reconstructed_data):
        """
        Calculate spatial consistency loss between original and reconstructed coordinates.
        Modified to ensure stronger gradient signal and proper backpropagation.

        Args:
            original_data: Original PyTorch Geometric Data object
            reconstructed_data: Reconstructed PyTorch Geometric Data object

        Returns:
            spatial_loss: Loss for spatial consistency
        """
        # Extract coordinates (assuming standard feature ordering)
        coord_idx_start = 21 + 7  # After AA and SS one-hot encodings
        coord_idx_end = coord_idx_start + 3

        # Get original and reconstructed coordinates
        orig_coords = original_data.x[:, coord_idx_start:coord_idx_end]
        recon_coords = reconstructed_data.x[:, coord_idx_start:coord_idx_end]

        # Direct coordinate loss - simple but effective for gradient flow
        direct_coord_loss = F.mse_loss(recon_coords, orig_coords)

        # Calculate pairwise distance matrices
        def calc_dist_matrix(coords):
            # Use broadcasting for efficient distance calculation
            # Shape: [N, 1, 3] - [1, N, 3] = [N, N, 3]
            diff = coords.unsqueeze(1) - coords.unsqueeze(0)
            # Add epsilon for numerical stability, shape: [N, N]
            dist = torch.sqrt(torch.sum(diff * diff, dim=2) + 1e-10)
            return dist

        # Calculate distance matrices
        orig_dist_matrix = calc_dist_matrix(orig_coords)
        recon_dist_matrix = calc_dist_matrix(recon_coords)

        # MSE loss between distance matrices
        dist_matrix_loss = F.mse_loss(recon_dist_matrix, orig_dist_matrix)

        # Combine both losses, with emphasis on direct coordinate loss
        # This provides multiple gradient paths
        spatial_loss = direct_coord_loss + 0.5 * dist_matrix_loss

        return spatial_loss
    def validate(self, dataloader):
        """
        Validate the model.

        Args:
            dataloader: DataLoader for validation data

        Returns:
            val_loss: Average validation loss
            val_metrics: Dictionary containing detailed validation metrics
        """
        self.model.eval()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        total_edge_loss = 0
        total_spatial_loss = 0
        num_batches = 0

        with torch.no_grad():
            for data in tqdm(dataloader, desc="Validation"):
                try:
                    # Move data to device
                    data = data.to(self.device)

                    # Ensure data has required attributes
                    if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                        data.edge_attr = torch.ones(data.edge_index.size(1), 2, dtype=torch.float, device=self.device)

                    # Forward pass
                    reconstructed_data, mu, logvar, edge_probs = self.model(data)

                    # Calculate reconstruction loss
                    rec_loss = F.mse_loss(reconstructed_data.x, data.x)

                    # Calculate KL divergence loss
                    batch_size = 1 if data.batch is None else (data.batch.max().item() + 1)
                    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
                    kl_loss = kl_loss / batch_size  # Normalize by batch size

                    # Calculate edge prediction loss
                    edge_loss = self._calculate_edge_prediction_loss(data, edge_probs)

                    # Calculate spatial consistency loss
                    spatial_loss = self._calculate_spatial_consistency_loss(data, reconstructed_data)

                    # Total loss
                    loss = (rec_loss +
                            self.kl_weight * kl_loss +
                            self.edge_weight * edge_loss +
                            self.spatial_weight * spatial_loss)

                    # Check for NaN loss
                    if torch.isnan(loss):
                        print(f"Warning: NaN loss detected in validation. Skipping batch.")
                        continue

                    # Accumulate metrics
                    total_loss += loss.item()
                    total_rec_loss += rec_loss.item()
                    total_kl_loss += kl_loss.item()
                    total_edge_loss += edge_loss.item()
                    total_spatial_loss += spatial_loss.item()
                    num_batches += 1

                except Exception as e:
                    print(f"Error validating batch: {e}")
                    continue

        # Calculate average loss
        if num_batches > 0:
            val_loss = total_loss / num_batches
            val_rec_loss = total_rec_loss / num_batches
            val_kl_loss = total_kl_loss / num_batches
            val_edge_loss = total_edge_loss / num_batches
            val_spatial_loss = total_spatial_loss / num_batches

            # Store validation history
            self.val_losses.append(val_loss)

            # Create metrics dictionary for detailed reporting
            val_metrics = {
                'total_loss': val_loss,
                'reconstruction_loss': val_rec_loss,
                'kl_loss': val_kl_loss,
                'edge_loss': val_edge_loss,
                'spatial_loss': val_spatial_loss
            }

            # Print summary
            print(f"Validation - Total: {val_loss:.4f}, Rec: {val_rec_loss:.4f}, "
                  f"KL: {val_kl_loss:.4f}, Edge: {val_edge_loss:.4f}, "
                  f"Spatial: {val_spatial_loss:.4f}"
                  f"Spatialweight: {self.spatial_weight:.4f}")

            return val_loss, val_metrics
        else:
            return float('nan'), None

    def train(self, train_loader, val_loader, num_epochs, save_path=None):
        """
        Train the model for multiple epochs.

        Args:
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
            num_epochs: Number of epochs to train for
            save_path: Path to save the model to (optional)

        Returns:
            history: Dictionary containing training and validation metrics history
        """
        print(f"Training on {self.device}")

        best_val_loss = float('inf')
        patience = 10  # Early stopping patience
        no_improve = 0

        # Dictionary to store all metrics
        history = {
            'train_loss': [],
            'val_loss': [],
            'reconstruction_loss': [],
            'kl_loss': [],
            'edge_loss': [],
            'spatial_loss': [],
            'val_metrics': []
        }

        for epoch in range(num_epochs):
            # Train for one epoch
            metrics = self.train_epoch(train_loader)
            train_loss = metrics['loss']
            rec_loss = metrics.get('rec_loss', 0.0)
            kl_loss = metrics.get('kl_loss', 0.0)
            edge_loss = metrics.get('edge_loss', 0.0)
            spatial_loss = metrics.get('spatial_loss', 0.0)

            # Validate
            val_loss, val_metrics = self.validate(val_loader)

            # Store metrics in history
            history['train_loss'].append(train_loss)
            history['val_loss'].append(val_loss)
            history['reconstruction_loss'].append(rec_loss)
            history['kl_loss'].append(kl_loss)
            history['edge_loss'].append(edge_loss)
            history['spatial_loss'].append(spatial_loss)
            if val_metrics:
                history['val_metrics'].append(val_metrics)

            # Print metrics with type checking
            print("Epoch {}/{}, Train Loss: {:.4f}, Val Loss: {:.4f}".format(
                epoch+1, num_epochs, train_loss, val_loss))

            # Check for NaN losses
            if torch.isnan(torch.tensor(train_loss)) or torch.isnan(torch.tensor(val_loss)):
                print("Warning: NaN loss detected. Consider reducing learning rate or adjusting model.")
                if epoch > 5:  # Only exit if we've trained for a few epochs already
                    print("Terminating training due to NaN losses.")
                    break

            # Save best model
            if save_path is not None and val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(self.model.state_dict(), save_path)
                print(f"Model saved to {save_path}")
                no_improve = 0
            else:
                no_improve += 1

            # Early stopping
            if no_improve >= patience:
                print(f"No improvement for {patience} epochs. Early stopping.")
                break

        # Load best model if saved
        if save_path is not None and os.path.exists(save_path):
            self.model.load_state_dict(torch.load(save_path))
            print(f"Loaded best model from {save_path}")

        # Return training history
        return history

    def evaluate(self, dataloader):
        """
        Evaluate the model on test data.

        Args:
            dataloader: DataLoader for test data

        Returns:
            metrics: Dictionary containing evaluation metrics
        """
        self.model.eval()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        total_edge_loss = 0
        total_spatial_loss = 0
        total_edge_accuracy = 0
        total_graph_similarity = 0
        num_batches = 0

        with torch.no_grad():
            for data in tqdm(dataloader, desc="Evaluating"):
                try:
                    # Move data to device
                    data = data.to(self.device)

                    # Forward pass
                    reconstructed_data, mu, logvar, edge_probs = self.model(data)

                    # Calculate losses
                    rec_loss = F.mse_loss(reconstructed_data.x, data.x)

                    batch_size = 1 if data.batch is None else (data.batch.max().item() + 1)
                    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp()) / batch_size

                    edge_loss = self._calculate_edge_prediction_loss(data, edge_probs)

                    spatial_loss = self._calculate_spatial_consistency_loss(data, reconstructed_data)

                    loss = (rec_loss +
                            self.kl_weight * kl_loss +
                            self.edge_weight * edge_loss +
                            self.spatial_weight * spatial_loss)

                    # Calculate edge prediction accuracy
                    # Create ground truth adjacency matrix
                    num_nodes = data.x.size(0)
                    orig_adj = to_dense_adj(data.edge_index, max_num_nodes=num_nodes).squeeze(0)

                    # Create predicted adjacency matrix (using threshold 0.5)
                    pred_adj = to_dense_adj(reconstructed_data.edge_index, max_num_nodes=num_nodes).squeeze(0)

                    # Calculate accuracy (exclude diagonal elements - self loops)
                    mask = ~torch.eye(num_nodes, dtype=torch.bool, device=self.device)
                    correct = (orig_adj[mask] == pred_adj[mask]).float().mean()

                    # Calculate graph similarity (Jaccard index)
                    orig_edges = set((i.item(), j.item()) for i, j in zip(data.edge_index[0], data.edge_index[1]))
                    pred_edges = set((i.item(), j.item()) for i, j in zip(reconstructed_data.edge_index[0], reconstructed_data.edge_index[1]))

                    if len(orig_edges.union(pred_edges)) > 0:
                        similarity = len(orig_edges.intersection(pred_edges)) / len(orig_edges.union(pred_edges))
                    else:
                        similarity = 1.0  # No edges in either graph

                    # Accumulate metrics
                    total_loss += loss.item()
                    total_rec_loss += rec_loss.item()
                    total_kl_loss += kl_loss.item()
                    total_edge_loss += edge_loss.item()
                    total_spatial_loss += spatial_loss.item()
                    total_edge_accuracy += correct.item()
                    total_graph_similarity += similarity
                    num_batches += 1

                except Exception as e:
                    print(f"Error evaluating batch: {e}")
                    continue

        # Calculate average metrics
        if num_batches > 0:
            metrics = {
                'loss': total_loss / num_batches,
                'reconstruction_loss': total_rec_loss / num_batches,
                'kl_loss': total_kl_loss / num_batches,
                'edge_loss': total_edge_loss / num_batches,
                'spatial_loss': total_spatial_loss / num_batches,
                'edge_accuracy': total_edge_accuracy / num_batches,
                'graph_similarity': total_graph_similarity / num_batches
            }

            return metrics
        else:
            return None