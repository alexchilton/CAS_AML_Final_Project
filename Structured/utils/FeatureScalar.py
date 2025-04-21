import numpy as np
import json
import os


class FeatureScaler:
    """
    Utility class for scaling and descaling protein graph features.
    Stores and manages scaling parameters for reproducible feature transformations.
    """

    def __init__(self, amino_acids=None, ss_types=None):
        """
        Initialize the feature scaler.

        Args:
            amino_acids: List of amino acid types (optional)
            ss_types: List of secondary structure types (optional)
        """
        # Feature type indices
        self.amino_acids = amino_acids or [
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL', 'UNK'
        ]
        self.ss_types = ss_types or ['H', 'E', 'B', 'G', 'T', 'S', '?']

        # Scaling parameters storage
        self.scaling_params = {
            'coordinates': {
                'mean': None,
                'std': None
            },
            'b_factor': {
                'mean': None,
                'std': None
            },
            'meiler_features': {
                'mean': None,
                'std': None
            }
        }

    def fit(self, features):
        """
        Compute scaling parameters from a set of features.

        Args:
            features: NumPy array of features to compute scaling from
        """
        # Coordinate indices
        coord_start = len(self.amino_acids) + len(self.ss_types)
        coord_end = coord_start + 3

        # B-factor index
        b_factor_idx = coord_end

        # Meiler features indices
        meiler_start = b_factor_idx + 1
        meiler_end = meiler_start + 7

        # Compute coordinate scaling
        coords = features[:, coord_start:coord_end]
        self.scaling_params['coordinates']['mean'] = np.mean(coords, axis=0)
        self.scaling_params['coordinates']['std'] = np.std(coords, axis=0)
        self.scaling_params['coordinates']['std'][self.scaling_params['coordinates']['std'] == 0] = 1.0

        # Compute B-factor scaling
        b_factors = features[:, b_factor_idx:b_factor_idx+1]
        self.scaling_params['b_factor']['mean'] = np.mean(b_factors)
        self.scaling_params['b_factor']['std'] = np.std(b_factors)
        if self.scaling_params['b_factor']['std'] == 0:
            self.scaling_params['b_factor']['std'] = 1.0

        # Compute Meiler features scaling
        meiler_features = features[:, meiler_start:meiler_end]
        self.scaling_params['meiler_features']['mean'] = np.mean(meiler_features, axis=0)
        self.scaling_params['meiler_features']['std'] = np.std(meiler_features, axis=0)
        self.scaling_params['meiler_features']['std'][self.scaling_params['meiler_features']['std'] == 0] = 1.0

    def transform(self, features):
        """
        Scale features using pre-computed parameters.

        Args:
            features: NumPy array of features to scale

        Returns:
            Scaled features
        """
        # Check if scaling parameters have been computed
        if any(param is None for group in self.scaling_params.values() for param in group.values()):
            raise ValueError("Scaling parameters not computed. Call fit() first.")

        # Create a copy to avoid modifying the original
        scaled_features = features.copy()

        # Coordinate indices
        coord_start = len(self.amino_acids) + len(self.ss_types)
        coord_end = coord_start + 3

        # B-factor index
        b_factor_idx = coord_end

        # Meiler features indices
        meiler_start = b_factor_idx + 1
        meiler_end = meiler_start + 7

        # Scale coordinates
        coords = scaled_features[:, coord_start:coord_end]
        scaled_features[:, coord_start:coord_end] = (
                (coords - self.scaling_params['coordinates']['mean']) /
                self.scaling_params['coordinates']['std']
        )

        # Scale B-factor
        b_factors = scaled_features[:, b_factor_idx:b_factor_idx+1]
        scaled_features[:, b_factor_idx:b_factor_idx+1] = (
                (b_factors - self.scaling_params['b_factor']['mean']) /
                self.scaling_params['b_factor']['std']
        )

        # Scale Meiler features
        meiler_features = scaled_features[:, meiler_start:meiler_end]
        scaled_features[:, meiler_start:meiler_end] = (
                (meiler_features - self.scaling_params['meiler_features']['mean']) /
                self.scaling_params['meiler_features']['std']
        )

        return scaled_features

    def inverse_transform(self, scaled_features):
        """
        Descale features back to original scale.

        Args:
            scaled_features: NumPy array of scaled features

        Returns:
            Descaled features
        """
        # Check if scaling parameters have been computed
        if any(param is None for group in self.scaling_params.values() for param in group.values()):
            raise ValueError("Scaling parameters not computed. Call fit() first.")

        # Create a copy to avoid modifying the original
        descaled_features = scaled_features.copy()

        # Coordinate indices
        coord_start = len(self.amino_acids) + len(self.ss_types)
        coord_end = coord_start + 3

        # B-factor index
        b_factor_idx = coord_end

        # Meiler features indices
        meiler_start = b_factor_idx + 1
        meiler_end = meiler_start + 7

        # Descale coordinates
        coords = descaled_features[:, coord_start:coord_end]
        descaled_features[:, coord_start:coord_end] = (
                coords * self.scaling_params['coordinates']['std'] +
                self.scaling_params['coordinates']['mean']
        )

        # Descale B-factor
        b_factors = descaled_features[:, b_factor_idx:b_factor_idx+1]
        descaled_features[:, b_factor_idx:b_factor_idx+1] = (
                b_factors * self.scaling_params['b_factor']['std'] +
                self.scaling_params['b_factor']['mean']
        )

        # Descale Meiler features
        meiler_features = descaled_features[:, meiler_start:meiler_end]
        descaled_features[:, meiler_start:meiler_end] = (
                meiler_features * self.scaling_params['meiler_features']['std'] +
                self.scaling_params['meiler_features']['mean']
        )

        # Convert continuous probabilities back to one-hot for amino acids and secondary structure
        # Amino acid features
        aa_features = descaled_features[:, :len(self.amino_acids)]
        aa_one_hot = np.zeros_like(aa_features)
        aa_one_hot[np.arange(len(aa_features)), aa_features.argmax(axis=1)] = 1
        descaled_features[:, :len(self.amino_acids)] = aa_one_hot

        # Secondary structure features
        ss_start = len(self.amino_acids)
        ss_end = ss_start + len(self.ss_types)
        ss_features = descaled_features[:, ss_start:ss_end]
        ss_one_hot = np.zeros_like(ss_features)
        ss_one_hot[np.arange(len(ss_features)), ss_features.argmax(axis=1)] = 1
        descaled_features[:, ss_start:ss_end] = ss_one_hot

        return descaled_features

    def save(self, filepath):
        """
        Save scaling parameters to a JSON file.

        Args:
            filepath: Path to save the scaling parameters
        """
        # Convert numpy arrays to lists for JSON serialization
        saveable_params = {}
        for group, params in self.scaling_params.items():
            saveable_params[group] = {
                key: value.tolist() if isinstance(value, np.ndarray) else value
                for key, value in params.items()
            }

        with open(filepath, 'w') as f:
            json.dump(saveable_params, f, indent=2)

    def load(self, filepath):
        """
        Load scaling parameters from a JSON file.

        Args:
            filepath: Path to load the scaling parameters from
        """
        with open(filepath, 'r') as f:
            loaded_params = json.load(f)

        # Convert lists back to numpy arrays
        for group, params in loaded_params.items():
            for key, value in params.items():
                self.scaling_params[group][key] = np.array(value)