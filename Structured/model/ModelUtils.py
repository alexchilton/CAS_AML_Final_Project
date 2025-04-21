import torch
import os
from typing import Dict, Any, Optional
import json

from GATVAE import GATVAE


class ModelUtils:
    """
    Utility class for saving and loading models and configurations.
    """

    @staticmethod
    def save_model(model: GATVAE, save_path: str, config: Optional[Dict[str, Any]] = None):
        """
        Save model weights and configuration.

        Args:
            model: GATVAE model to save
            save_path: Path to save model weights
            config: Optional configuration dictionary to save with the model
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(save_path), exist_ok=True)

        # Save model weights
        torch.save(model.state_dict(), save_path)

        # Save config if provided
        if config:
            config_path = save_path.replace('.pth', '_config.json')
            with open(config_path, 'w') as f:
                json.dump(config, f, indent=2)

        print(f"Model saved to {save_path}")

    @staticmethod
    def load_model(load_path: str, input_dim: int, config_path: Optional[str] = None) -> GATVAE:
        """
        Load model weights and configuration.

        Args:
            load_path: Path to load model weights from
            input_dim: Input dimension for model initialization
            config_path: Optional path to load configuration from

        Returns:
            Loaded GATVAE model
        """
        # Load config if provided
        if config_path is None:
            config_path = load_path.replace('.pth', '_config.json')

        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                config = json.load(f)

            # Initialize model with loaded configuration
            model = GATVAE(
                input_dim=input_dim,
                hidden_dim=config.get('hidden_dim', 64),
                latent_dim=config.get('latent_dim', 32),
                num_heads=config.get('num_heads', 4),
                dropout=config.get('dropout', 0.1)
            )
        else:
            # Initialize model with default parameters
            model = GATVAE(
                input_dim=input_dim,
                hidden_dim=64,
                latent_dim=32,
                num_heads=4,
                dropout=0.1
            )

        # Load model weights
        model.load_state_dict(torch.load(load_path))
        model.eval()

        print(f"Model loaded from {load_path}")

        return model

    @staticmethod
    def export_latent_vectors(latent_vectors, ids, export_path: str):
        """
        Export latent vectors and their corresponding IDs to a file.

        Args:
            latent_vectors: List of latent vectors
            ids: List of corresponding IDs
            export_path: Path to save the latent vectors to
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(export_path), exist_ok=True)

        # Convert latent vectors to numpy
        import numpy as np
        if isinstance(latent_vectors, list):
            latent_vectors = np.vstack(latent_vectors)

        # Save to file
        np.savez(
            export_path,
            latent_vectors=latent_vectors,
            ids=ids
        )

        print(f"Latent vectors exported to {export_path}")

    @staticmethod
    def import_latent_vectors(import_path: str):
        """
        Import latent vectors and their corresponding IDs from a file.

        Args:
            import_path: Path to load the latent vectors from

        Returns:
            Tuple of (latent_vectors, ids)
        """
        import numpy as np

        # Load from file
        data = np.load(import_path)
        latent_vectors = data['latent_vectors']
        ids = data['ids']

        print(f"Latent vectors imported from {import_path}")

        return latent_vectors, ids