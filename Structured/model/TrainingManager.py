import torch
import torch.nn.functional as F
from torch.optim import Adam
from torch.utils.data import DataLoader
from torch_geometric.data import Data, Batch
import numpy as np
from tqdm import tqdm
import os

from GATVAE import GATVAE


class TrainingManager:
    """
    Manager for training the GAT-based VAE model.
    Handles training loop, loss calculation, and evaluation.
    """

    def __init__(self, model, learning_rate=0.001, kl_weight=0.01):
        """
        Initialize the training manager.

        Args:
            model: GATVAE model to train
            learning_rate: Learning rate for optimization
            kl_weight: Weight for the KL divergence loss term
        """
        self.model = model
        self.optimizer = Adam(model.parameters(), lr=learning_rate)
        self.kl_weight = kl_weight
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Move model to device
        self.model.to(self.device)

        # Training history
        self.train_losses = []
        self.val_losses = []
        self.reconstruction_losses = []
        self.kl_losses = []

    def train_epoch(self, dataloader):
        """
        Train for one epoch.

        Args:
            dataloader: DataLoader for training data

        Returns:
            epoch_loss: Average loss for the epoch
        """
        self.model.train()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        num_batches = 0

        for data in tqdm(dataloader, desc="Training"):
            try:
                # Move data to device
                data = data.to(self.device)

                # Ensure data has required attributes
                if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                    data.edge_attr = torch.ones(data.edge_index.size(1), 1, dtype=torch.float, device=self.device)

                # Zero gradients
                self.optimizer.zero_grad()

                # Forward pass
                x_reconstructed, mu, logvar = self.model(data)

                # Calculate reconstruction loss
                # Ensure proper scaling for batch size
                rec_loss = F.mse_loss(x_reconstructed, data.x)

                # Calculate KL divergence loss with proper scaling
                # Scale KL loss by the number of graphs in the batch
                batch_size = 1 if data.batch is None else (data.batch.max().item() + 1)
                kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
                kl_loss = kl_loss / batch_size  # Normalize by batch size

                # Total loss
                loss = rec_loss + self.kl_weight * kl_loss

                # Check for NaN loss
                if torch.isnan(loss):
                    print(f"Warning: NaN loss detected. Skipping batch.")
                    print(f"  Reconstruction loss: {rec_loss.item()}")
                    print(f"  KL loss: {kl_loss.item()}")
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
        else:
            epoch_loss = float('nan')
            epoch_rec_loss = float('nan')
            epoch_kl_loss = float('nan')

        # Store training history
        self.train_losses.append(epoch_loss)
        self.reconstruction_losses.append(epoch_rec_loss)
        self.kl_losses.append(epoch_kl_loss)

        return epoch_loss, epoch_rec_loss, epoch_kl_loss
    def validate(self, dataloader):
        """
        Validate the model.

        Args:
            dataloader: DataLoader for validation data

        Returns:
            val_loss: Average validation loss
        """
        self.model.eval()
        total_loss = 0
        total_rec_loss = 0
        total_kl_loss = 0
        num_batches = 0

        with torch.no_grad():
            for data in tqdm(dataloader, desc="Validation"):
                try:
                    # Move data to device
                    data = data.to(self.device)

                    # Ensure data has required attributes
                    if not hasattr(data, 'edge_attr') or data.edge_attr is None:
                        data.edge_attr = torch.ones(data.edge_index.size(1), 1, dtype=torch.float, device=self.device)

                    # Forward pass
                    x_reconstructed, mu, logvar = self.model(data)

                    # Calculate reconstruction loss
                    rec_loss = F.mse_loss(x_reconstructed, data.x)

                    # Calculate KL divergence loss with proper scaling
                    batch_size = 1 if data.batch is None else (data.batch.max().item() + 1)
                    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
                    kl_loss = kl_loss / batch_size  # Normalize by batch size

                    # Total loss
                    loss = rec_loss + self.kl_weight * kl_loss

                    # Check for NaN loss
                    if torch.isnan(loss):
                        print(f"Warning: NaN loss detected in validation. Skipping batch.")
                        continue

                    # Accumulate metrics
                    total_loss += loss.item()
                    total_rec_loss += rec_loss.item()
                    total_kl_loss += kl_loss.item()
                    num_batches += 1

                except Exception as e:
                    print(f"Error validating batch: {e}")
                    continue

        # Calculate average loss
        if num_batches > 0:
            val_loss = total_loss / num_batches
            val_rec_loss = total_rec_loss / num_batches
            val_kl_loss = total_kl_loss / num_batches

            print(f"Validation - Rec Loss: {val_rec_loss:.4f}, KL Loss: {val_kl_loss:.4f}")
        else:
            val_loss = float('nan')

        # Store validation history
        self.val_losses.append(val_loss)

        return val_loss

    def train(self, train_loader, val_loader, num_epochs, save_path=None):
        """
        Train the model for multiple epochs.

        Args:
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
            num_epochs: Number of epochs to train for
            save_path: Path to save the model to (optional)

        Returns:
            Training history
        """
        print(f"Training on {self.device}")

        best_val_loss = float('inf')
        patience = 5  # Early stopping patience
        no_improve = 0

        for epoch in range(num_epochs):
            # Train for one epoch
            train_loss, rec_loss, kl_loss = self.train_epoch(train_loader)

            # Validate
            val_loss = self.validate(val_loader)

            # Print metrics
            print(f"Epoch {epoch+1}/{num_epochs}, "
                  f"Train Loss: {train_loss:.4f}, "
                  f"Rec Loss: {rec_loss:.4f}, "
                  f"KL Loss: {kl_loss:.4f}, "
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
        return {
            'train_losses': self.train_losses,
            'val_losses': self.val_losses,
            'reconstruction_losses': self.reconstruction_losses,
            'kl_losses': self.kl_losses
        }
    def evaluate(self, dataloader):
        """
        Evaluate the model on test data.

        Args:
            dataloader: DataLoader for test data

        Returns:
            test_loss: Average test loss
        """
        self.model.eval()
        total_loss = 0

        with torch.no_grad():
            for data in tqdm(dataloader, desc="Testing"):
                # Move data to device
                data = data.to(self.device)

                # Forward pass
                x_reconstructed, mu, logvar = self.model(data)

                # Calculate loss
                rec_loss = F.mse_loss(x_reconstructed, data.x)
                kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
                loss = rec_loss + self.kl_weight * kl_loss

                # Accumulate metrics
                total_loss += loss.item()

        # Calculate average loss
        test_loss = total_loss / len(dataloader)

        return test_loss