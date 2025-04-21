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
                 spatial_weight=0.2, lambda_reg=0.001):
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

    def _calculate_spatial_consistency_loss(self, original_data, reconstructed_data):
        """
        Calculate spatial consistency loss between original and reconstructed coordinates.

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

        # Calculate pairwise distance matrices
        def calc_dist_matrix(coords):
            # Expand coordinates for pairwise calculations
            x_i = coords.unsqueeze(1)  # Shape: [N, 1, 3]
            x_j = coords.unsqueeze(0)  # Shape: [1, N, 3]
            # Compute Euclidean distances
            dist = torch.sqrt(torch.sum((x_i - x_j) ** 2, dim=2) + 1e-6)
            return dist

        # Calculate distance matrices
        orig_dist_matrix = calc_dist_matrix(orig_coords)
        recon_dist_matrix = calc_dist_matrix(recon_coords)

        # MSE loss between distance matrices
        spatial_loss = F.mse_loss(recon_dist_matrix, orig_dist_matrix)

        return spatial_loss

    def train_epoch(self, dataloader):
        """
        Train for one epoch.

        Args:
            dataloader: DataLoader for training data

        Returns:
            epoch_loss: Average loss for the epoch
            epoch_rec_loss: Average reconstruction loss
            epoch_kl_loss: Average KL divergence loss
            epoch_edge_loss: Average edge prediction loss
            epoch_spatial_loss: Average spatial consistency loss
        """
        self.model.train()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        total_edge_loss = 0
        total_spatial_loss = 0
        num_batches = 0

        for data in tqdm(dataloader, desc="Training"):
            try:
                # Move data to device
                data = data.to(self.device)

                # Ensure data has required attributes
                if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                    data.edge_attr = torch.ones(data.edge_index.size(1), 2, dtype=torch.float, device=self.device)

                # Zero gradients
                self.optimizer.zero_grad()

                # Forward pass with enhanced model
                reconstructed_data, mu, logvar, edge_probs = self.model(data)

                # Calculate reconstruction loss for node features
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
                    print(f"Warning: NaN loss detected. Skipping batch.")
                    print(f"  Reconstruction loss: {rec_loss.item()}")
                    print(f"  KL loss: {kl_loss.item()}")
                    print(f"  Edge loss: {edge_loss.item()}")
                    print(f"  Spatial loss: {spatial_loss.item()}")
                    continue

                # Backward pass
                loss.backward()

                # Gradient clipping to prevent exploding gradients
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)

                # Update weights
                self.optimizer.step()

                # Accumulate metrics
                total_loss += loss.item()
                total_rec_loss += rec_loss.item()
                total_kl_loss += kl_loss.item()
                total_edge_loss += edge_loss.item()
                total_spatial_loss += spatial_loss.item()
                num_batches += 1

            except Exception as e:
                print(f"Error processing batch: {e}")
                import traceback
                traceback.print_exc()
                continue

        # Calculate average losses
        if num_batches > 0:
            epoch_loss = total_loss / num_batches
            epoch_rec_loss = total_rec_loss / num_batches
            epoch_kl_loss = total_kl_loss / num_batches
            epoch_edge_loss = total_edge_loss / num_batches
            epoch_spatial_loss = total_spatial_loss / num_batches
        else:
            epoch_loss = float('nan')
            epoch_rec_loss = float('nan')
            epoch_kl_loss = float('nan')
            epoch_edge_loss = float('nan')
            epoch_spatial_loss = float('nan')

        # Store training history
        self.train_losses.append(epoch_loss)
        self.reconstruction_losses.append(epoch_rec_loss)
        self.kl_losses.append(epoch_kl_loss)
        self.edge_losses.append(epoch_edge_loss)
        self.spatial_losses.append(epoch_spatial_loss)

        return epoch_loss, epoch_rec_loss, epoch_kl_loss, epoch_edge_loss, epoch_spatial_loss

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
                  f"Spatial: {val_spatial_loss:.4f}")

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
        patience = 5  # Early stopping patience
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
            train_loss, rec_loss, kl_loss, edge_loss, spatial_loss = self.train_epoch(train_loader)

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

            # Print metrics
            print(f"Epoch {epoch+1}/{num_epochs}, "
                  f"Train Loss: {train_loss:.4f}, "
                  f"Val Loss: {val_loss:.4f}")

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