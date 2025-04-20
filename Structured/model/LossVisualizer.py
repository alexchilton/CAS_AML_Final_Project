import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Any


class LossVisualizer:
    """
    A simple class to visualize training and validation losses.
    """

    def __init__(self, save_dir='./plots'):
        """
        Initialize the loss visualizer.

        Args:
            save_dir: Directory to save plots
        """
        self.save_dir = save_dir

        # Create directory if it doesn't exist
        import os
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

    def plot_training_history(self, history: Dict[str, List[float]], show=True, save=True):
        """
        Plot training and validation losses.

        Args:
            history: Dictionary containing loss histories
            show: Whether to display the plot
            save: Whether to save the plot to disk
        """
        epochs = range(1, len(history['train_losses']) + 1)

        plt.figure(figsize=(10, 6))
        plt.plot(epochs, history['train_losses'], 'b-', label='Training Loss')
        plt.plot(epochs, history['val_losses'], 'r-', label='Validation Loss')
        plt.title('Training and Validation Loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.grid(True)

        if save:
            plt.savefig(f"{self.save_dir}/train_val_loss.png", dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

    def plot_component_losses(self, history: Dict[str, List[float]], show=True, save=True):
        """
        Plot reconstruction and KL divergence losses.

        Args:
            history: Dictionary containing loss histories
            show: Whether to display the plot
            save: Whether to save the plot to disk
        """
        epochs = range(1, len(history['reconstruction_losses']) + 1)

        plt.figure(figsize=(10, 6))
        plt.plot(epochs, history['reconstruction_losses'], 'g-', label='Reconstruction Loss')
        plt.plot(epochs, history['kl_losses'], 'p-', label='KL Divergence Loss')
        plt.title('Component Losses')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.legend()
        plt.grid(True)

        if save:
            plt.savefig(f"{self.save_dir}/component_losses.png", dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

    def plot_loss_ratio(self, history: Dict[str, List[float]], show=True, save=True):
        """
        Plot the ratio of reconstruction loss to KL divergence loss.

        Args:
            history: Dictionary containing loss histories
            show: Whether to display the plot
            save: Whether to save the plot to disk
        """
        epochs = range(1, len(history['reconstruction_losses']) + 1)

        # Calculate ratio
        ratios = [r / k if k != 0 else 0 for r, k in
                  zip(history['reconstruction_losses'], history['kl_losses'])]

        plt.figure(figsize=(10, 6))
        plt.plot(epochs, ratios, 'c-', label='Rec Loss / KL Loss')
        plt.title('Reconstruction to KL Loss Ratio')
        plt.xlabel('Epochs')
        plt.ylabel('Ratio')
        plt.legend()
        plt.grid(True)

        if save:
            plt.savefig(f"{self.save_dir}/loss_ratio.png", dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

    def plot_all(self, history: Dict[str, List[float]], show=True, save=True):
        """
        Create all available plots.

        Args:
            history: Dictionary containing loss histories
            show: Whether to display the plots
            save: Whether to save the plots to disk
        """
        self.plot_training_history(history, show, save)
        self.plot_component_losses(history, show, save)
        self.plot_loss_ratio(history, show, save)