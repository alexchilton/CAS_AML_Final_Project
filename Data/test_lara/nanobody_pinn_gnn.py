import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATv2Conv, global_mean_pool, global_add_pool
from torch_geometric.data import DataLoader, Dataset
import numpy as np
import json
from tqdm import tqdm

class NanobodyDataset(Dataset):
    """PyTorch Geometric Dataset for Nanobody graphs."""
    
    def __init__(self, file_paths):
        """
        Initialize the dataset.
        
        Parameters:
        -----------
        file_paths : list
            List of paths to .pt files containing PyTorch Geometric graphs
        """
        super(NanobodyDataset, self).__init__()
        self.file_paths = file_paths
    
    def len(self):
        return len(self.file_paths)
    
    def get(self, idx):
        """Load a graph from disk."""
        return torch.load(self.file_paths[idx])


class PhysicsInformedLayer(nn.Module):
    """
    A physics-informed layer that incorporates physical constraints
    into the neural network.
    """
    
    def __init__(self, in_channels, out_channels):
        """
        Initialize the physics-informed layer.
        
        Parameters:
        -----------
        in_channels : int
            Number of input channels
        out_channels : int
            Number of output channels
        """
        super(PhysicsInformedLayer, self).__init__()
        
        # Regular neural network part
        self.nn_transform = nn.Sequential(
            nn.Linear(in_channels, 2 * out_channels),
            nn.ReLU(),
            nn.Linear(2 * out_channels, out_channels)
        )
        
        # Physics-based transformation
        # For example, we might apply Henderson-Hasselbalch equation for pH
        # or Debye-Hückel theory for ionic strength
        self.charge_weight = nn.Parameter(torch.randn(1))
        self.hydrophobic_weight = nn.Parameter(torch.randn(1))
    
    def forward(self, x, charge, hydrophobic):
        """
        Forward pass with both NN and physics-based components.
        
        Parameters:
        -----------
        x : torch.Tensor
            Input features
        charge : torch.Tensor
            Charge information
        hydrophobic : torch.Tensor
            Hydrophobicity information
            
        Returns:
        --------
        torch.Tensor
            Transformed features
        """
        # Neural network transformation
        nn_output = self.nn_transform(x)
        
        # Physics-based transformation (simplified example)
        # In a real implementation, you would use more complex physical equations
        physics_term = self.charge_weight * charge + self.hydrophobic_weight * hydrophobic
        physics_output = physics_term.unsqueeze(-1).expand_as(nn_output)
        
        # Combine both transformations
        return nn_output + 0.1 * physics_output


class EnhancedPhysicsInformedLayer(nn.Module):
    """
    An enhanced physics-informed layer that incorporates more detailed
    physical constraints into the neural network.
    """
    
    def __init__(self, in_channels, out_channels):
        """
        Initialize the enhanced physics-informed layer.
        
        Parameters:
        -----------
        in_channels : int
            Number of input channels
        out_channels : int
            Number of output channels
        """
        super(EnhancedPhysicsInformedLayer, self).__init__()
        
        # Regular neural network part with residual connection
        self.nn_transform = nn.Sequential(
            nn.Linear(in_channels, 2 * out_channels),
            nn.LayerNorm(2 * out_channels),
            nn.SiLU(),  # SiLU (Swish) activation for better gradient flow
            nn.Linear(2 * out_channels, out_channels)
        )
        
        # Residual connection
        self.residual = nn.Linear(in_channels, out_channels) if in_channels != out_channels else nn.Identity()
        
        # Physics-based transformations
        # Henderson-Hasselbalch equation components
        self.pka_transform = nn.Linear(1, out_channels)
        self.charge_transform = nn.Linear(1, out_channels)
        
        # Debye-Hückel theory components
        self.hydrophobic_transform = nn.Linear(1, out_channels)
        self.ionic_transform = nn.Linear(1, out_channels)
        
        # Dynamic weighting for physics components
        self.physics_weights = nn.Parameter(torch.ones(4))
        self.physics_scale = nn.Parameter(torch.tensor(0.1))
        
        # Gating mechanism to control physics influence
        self.gate = nn.Sequential(
            nn.Linear(in_channels, out_channels),
            nn.Sigmoid()
        )
    
    def forward(self, x, charge, hydrophobic, pka=None, ionic_radius=None):
        """
        Forward pass with enhanced physics-based components.
        
        Parameters:
        -----------
        x : torch.Tensor
            Input features
        charge : torch.Tensor
            Charge information
        hydrophobic : torch.Tensor
            Hydrophobicity information
        pka : torch.Tensor, optional
            pKa values if available
        ionic_radius : torch.Tensor, optional
            Ionic radius if available
            
        Returns:
        --------
        torch.Tensor
            Transformed features
        """
        # Residual neural network pathway
        nn_output = self.nn_transform(x) + self.residual(x)
        
        # Prepare physics features
        charge_expanded = charge.unsqueeze(-1)
        hydrophobic_expanded = hydrophobic.unsqueeze(-1)
        
        # Default values for optional parameters
        if pka is None:
            pka = torch.ones_like(charge) * 7.0  # Default neutral pKa
        if ionic_radius is None:
            ionic_radius = torch.ones_like(charge) * 0.5  # Default medium ionic radius
            
        pka_expanded = pka.unsqueeze(-1)
        ionic_radius_expanded = ionic_radius.unsqueeze(-1)
        
        # Apply physics transformations
        ph_component = self.pka_transform(pka_expanded)
        charge_component = self.charge_transform(charge_expanded)
        hydrophobic_component = self.hydrophobic_transform(hydrophobic_expanded)
        ionic_component = self.ionic_transform(ionic_radius_expanded)
        
        # Combine physics components with learned weights
        physics_components = [
            ph_component * F.softplus(self.physics_weights[0]),
            charge_component * F.softplus(self.physics_weights[1]),
            hydrophobic_component * F.softplus(self.physics_weights[2]),
            ionic_component * F.softplus(self.physics_weights[3])
        ]
        
        physics_output = sum(physics_components)
        
        # Scale physics contribution
        physics_contribution = physics_output * F.softplus(self.physics_scale)
        
        # Calculate gating factor - determines how much physics to incorporate
        gate_value = self.gate(x)
        
        # Combine using gating mechanism
        output = nn_output + gate_value * physics_contribution
        
        return output


class NanobodyPINNGNN(nn.Module):
    """
    Physics-Informed Graph Neural Network for predicting optimal
    pH and ionic strength for nanobodies.
    """
    
    def __init__(self, node_features, edge_features, hidden_dim=64, num_layers=3, 
                 dropout=0.2, use_gat=True):
        """
        Initialize the PINN GNN model.
        
        Parameters:
        -----------
        node_features : int
            Number of node features
        edge_features : int
            Number of edge features
        hidden_dim : int
            Dimension of hidden layers
        num_layers : int
            Number of GNN layers
        dropout : float
            Dropout probability
        use_gat : bool
            Whether to use GAT instead of GCN
        """
        super(NanobodyPINNGNN, self).__init__()
        
        self.node_features = node_features
        self.edge_features = edge_features
        self.hidden_dim = hidden_dim
        self.use_gat = use_gat
        
        # Edge feature transformation
        self.edge_encoder = nn.Linear(edge_features, hidden_dim)
        
        # Graph convolution layers
        self.convs = nn.ModuleList()
        
        # Initialize the first layer with the number of node features
        if use_gat:
            self.convs.append(GATv2Conv(node_features, hidden_dim, edge_dim=hidden_dim))
        else:
            self.convs.append(GCNConv(node_features, hidden_dim))
        
        # Add additional layers
        for _ in range(num_layers - 1):
            if use_gat:
                self.convs.append(GATv2Conv(hidden_dim, hidden_dim, edge_dim=hidden_dim))
            else:
                self.convs.append(GCNConv(hidden_dim, hidden_dim))
        
        # Physics-informed layer
        self.physics_layer = PhysicsInformedLayer(hidden_dim, hidden_dim)
        
        # Pooling
        self.pool = global_mean_pool
        
        # Additional physics-informed components for pH prediction
        # Henderson-Hasselbalch inspired module
        self.ph_physics = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)  # pH output
        )
        
        # Debye-Hückel inspired module for ionic strength
        self.ionic_physics = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1)  # Ionic strength output
        )
        
        # Output layers (combines both neural network and physics components)
        self.fc1 = nn.Linear(hidden_dim, hidden_dim)
        self.fc2_ph = nn.Linear(hidden_dim, 1)  # pH output from pure NN
        self.fc2_ionic = nn.Linear(hidden_dim, 1)  # Ionic strength from pure NN
        
        # Physics-based constraint modules
        self.ph_constraint = nn.Sequential(
            nn.Linear(1, hidden_dim // 4),
            nn.Tanh(),
            nn.Linear(hidden_dim // 4, 1)
        )
        
        self.ionic_constraint = nn.Sequential(
            nn.Linear(1, hidden_dim // 4),
            nn.Tanh(),
            nn.Linear(hidden_dim // 4, 1)
        )
        
        # Blending parameters (learnable)
        self.alpha_ph = nn.Parameter(torch.tensor([0.5]))
        self.alpha_ionic = nn.Parameter(torch.tensor([0.5]))
        
        # Dropout
        self.dropout = nn.Dropout(dropout)
        
        # Batch normalization
        self.batch_norm = nn.BatchNorm1d(hidden_dim)
    
    def extract_physics_features(self, x):
        """
        Extract physics-relevant features from node features.
        
        Parameters:
        -----------
        x : torch.Tensor
            Node features
            
        Returns:
        --------
        tuple
            (charge, hydrophobicity)
        """
        # Extract features based on your node feature structure
        # This depends on how you organized your features in NanobodyGraphBuilder
        
        # Assuming:
        # - Mean charge is the second-to-last feature
        # - Hydrophobicity is the third-to-last feature
        charge = x[:, -2]
        hydrophobic = x[:, -3]
        
        return charge, hydrophobic
    
    def forward(self, data):
        """
        Forward pass through the model.
        
        Parameters:
        -----------
        data : torch_geometric.data.Data
            Input graph data
            
        Returns:
        --------
        torch.Tensor
            Predicted [pH, ionic_strength]
        """
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        
        # Extract physics features
        charge, hydrophobic = self.extract_physics_features(x)
        
        # Transform edge features
        edge_attr = self.edge_encoder(edge_attr)
        
        # Apply GNN layers
        for i, conv in enumerate(self.convs):
            if self.use_gat:
                x = conv(x, edge_index, edge_attr)
            else:
                x = conv(x, edge_index)
            
            x = F.relu(x)
            x = self.dropout(x)
        
        # Apply physics-informed layer
        x = self.physics_layer(x, charge, hydrophobic)
        
        # Global pooling
        x_pooled = self.pool(x, batch)
        
        # Apply batch normalization to pooled features
        x_pooled = self.batch_norm(x_pooled)
        
        # Compute shared representation
        shared_features = F.relu(self.fc1(x_pooled))
        shared_features = self.dropout(shared_features)
        
        # Pure neural network pathway
        ph_nn = self.fc2_ph(shared_features)
        ionic_nn = self.fc2_ionic(shared_features)
        
        # Physics-informed pathway
        ph_phys = self.ph_physics(x_pooled)
        ionic_phys = self.ionic_physics(x_pooled)
        
        # Calculate aggregate charge and hydrophobicity for the whole protein
        agg_charge = global_add_pool(charge.unsqueeze(1), batch)
        agg_hydrophobic = global_add_pool(hydrophobic.unsqueeze(1), batch)
        
        # Apply physics constraints
        # For pH, use Henderson-Hasselbalch inspired constraint
        # Log(acid/base) ≈ pH - pKa
        charge_ratio = torch.sigmoid(agg_charge) * 10.0  # Scale to reasonable pH range
        ph_constraint = self.ph_constraint(charge_ratio)
        
        # For ionic strength, use Debye-Hückel inspired constraint
        # Higher overall charge requires higher ionic strength
        ionic_constraint = self.ionic_constraint(torch.abs(agg_charge))
        
        # Blend neural network and physics-based predictions
        # The alpha parameters are learnable and determine the weight of each component
        ph_pred = self.alpha_ph * ph_nn + (1 - self.alpha_ph) * ph_phys + 0.1 * ph_constraint
        ionic_pred = self.alpha_ionic * ionic_nn + (1 - self.alpha_ionic) * ionic_phys + 0.1 * ionic_constraint
        
        # Ensure predictions are in valid ranges
        # pH typically ranges from 0 to 14
        ph_pred = 7.0 + torch.tanh(ph_pred) * 7.0  # Center around 7 with range 0-14
        
        # Ionic strength is non-negative (typically 0 to 1 mol/L for most buffers)
        ionic_pred = F.softplus(ionic_pred)  # Ensures positive values
        
        # Combine predictions
        predictions = torch.cat([ph_pred, ionic_pred], dim=1)
        
        return predictions


class NanobodyPINNTrainer:
    """Trainer for the Nanobody PINN GNN model."""
    
    def __init__(self, model, optimizer, scheduler=None, device=None):
        """
        Initialize the trainer.
        
        Parameters:
        -----------
        model : NanobodyPINNGNN
            The model to train
        optimizer : torch.optim.Optimizer
            The optimizer to use
        scheduler : torch.optim.lr_scheduler._LRScheduler, optional
            Learning rate scheduler
        device : torch.device, optional
            Device to use for training
        """
        self.model = model
        self.optimizer = optimizer
        self.scheduler = scheduler
        self.device = device if device is not None else torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        
        # Training history
        self.history = {
            'train_loss': [],
            'val_loss': [],
            'test_loss': [],
            'best_val_loss': float('inf')
        }
    
    def train_epoch(self, train_loader):
        """
        Train for one epoch.
        
        Parameters:
        -----------
        train_loader : torch_geometric.data.DataLoader
            DataLoader for training data
            
        Returns:
        --------
        float
            Average epoch loss
        """
        self.model.train()
        total_loss = 0
        
        for data in train_loader:
            # Move data to device
            data = data.to(self.device)
            
            # Zero gradients
            self.optimizer.zero_grad()
            
            # Forward pass
            out = self.model(data)
            
            # Calculate loss (MSE for pH and ionic strength)
            loss = F.mse_loss(out, data.y)
            
            # Backward pass
            loss.backward()
            
            # Update parameters
            self.optimizer.step()
            
            # Accumulate loss
            total_loss += loss.item() * data.num_graphs
        
        # Calculate average loss
        epoch_loss = total_loss / len(train_loader.dataset)
        
        # Update learning rate if scheduler exists
        if self.scheduler is not None:
            self.scheduler.step()
        
        return epoch_loss
    
    def evaluate(self, loader):
        """
        Evaluate the model on a dataset.
        
        Parameters:
        -----------
        loader : torch_geometric.data.DataLoader
            DataLoader for evaluation data
            
        Returns:
        --------
        tuple
            (average loss, predictions, targets)
        """
        self.model.eval()
        total_loss = 0
        all_preds = []
        all_targets = []
        
        with torch.no_grad():
            for data in loader:
                # Move data to device
                data = data.to(self.device)
                
                # Forward pass
                out = self.model(data)
                
                # Calculate loss
                loss = F.mse_loss(out, data.y)
                
                # Accumulate loss
                total_loss += loss.item() * data.num_graphs
                
                # Store predictions and targets
                all_preds.append(out.cpu())
                all_targets.append(data.y.cpu())
        
        # Calculate average loss
        avg_loss = total_loss / len(loader.dataset)
        
        # Concatenate predictions and targets
        all_preds = torch.cat(all_preds, dim=0)
        all_targets = torch.cat(all_targets, dim=0)
        
        return avg_loss, all_preds, all_targets
    
    def train(self, train_loader, val_loader, test_loader=None, 
              num_epochs=100, patience=10, checkpoint_dir='checkpoints'):
        """
        Train the model.
        
        Parameters:
        -----------
        train_loader : torch_geometric.data.DataLoader
            DataLoader for training data
        val_loader : torch_geometric.data.DataLoader
            DataLoader for validation data
        test_loader : torch_geometric.data.DataLoader, optional
            DataLoader for test data
        num_epochs : int
            Number of epochs to train for
        patience : int
            Early stopping patience
        checkpoint_dir : str
            Directory to save checkpoints
            
        Returns:
        --------
        dict
            Training history
        """
        os.makedirs(checkpoint_dir, exist_ok=True)
        
        best_val_loss = float('inf')
        patience_counter = 0
        
        for epoch in range(1, num_epochs + 1):
            # Train for one epoch
            train_loss = self.train_epoch(train_loader)
            
            # Evaluate on validation set
            val_loss, val_preds, val_targets = self.evaluate(val_loader)
            
            # Print progress
            print(f'Epoch {epoch}/{num_epochs}, Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}')
            
            # Update history
            self.history['train_loss'].append(train_loss)
            self.history['val_loss'].append(val_loss)
            
            # Check for improvement
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                
                # Save best model
                self.history['best_val_loss'] = best_val_loss
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': self.model.state_dict(),
                    'optimizer_state_dict': self.optimizer.state_dict(),
                    'val_loss': val_loss,
                    'train_loss': train_loss
                }, os.path.join(checkpoint_dir, 'best_model.pth'))
                
                print(f'Best model saved with validation loss: {best_val_loss:.4f}')
            else:
                patience_counter += 1
            
            # Early stopping
            if patience_counter >= patience:
                print(f'Early stopping after {epoch} epochs')
                break
        
        # Evaluate on test set if provided
        if test_loader is not None:
            # Load best model
            checkpoint = torch.load(os.path.join(checkpoint_dir, 'best_model.pth'))
            self.model.load_state_dict(checkpoint['model_state_dict'])
            
            # Evaluate
            test_loss, test_preds, test_targets = self.evaluate(test_loader)
            self.history['test_loss'] = test_loss
            
            print(f'Test Loss: {test_loss:.4f}')
            
            # Calculate metrics
            test_metrics = {
                'rmse': torch.sqrt(F.mse_loss(test_preds, test_targets)).item(),
                'mae': F.l1_loss(test_preds, test_targets).item(),
                'r2': self._r2_score(test_preds, test_targets)
            }
            
            # Separate metrics for pH and ionic strength
            test_metrics['ph_rmse'] = torch.sqrt(F.mse_loss(test_preds[:, 0], test_targets[:, 0])).item()
            test_metrics['ionic_rmse'] = torch.sqrt(F.mse_loss(test_preds[:, 1], test_targets[:, 1])).item()
            
            self.history['test_metrics'] = test_metrics
            
            print(f'Test RMSE: {test_metrics["rmse"]:.4f}, R²: {test_metrics["r2"]:.4f}')
            print(f'pH RMSE: {test_metrics["ph_rmse"]:.4f}, Ionic Strength RMSE: {test_metrics["ionic_rmse"]:.4f}')
        
        return self.history
    
    def _r2_score(self, pred, target):
        """Calculate R² score."""
        total_variance = torch.var(target, unbiased=False) * target.shape[0]
        explained_variance = torch.sum((pred - target) ** 2)
        r2 = 1 - (explained_variance / total_variance)
        return r2.item()
    
    def predict(self, loader):
        """
        Make predictions with the model.
        
        Parameters:
        -----------
        loader : torch_geometric.data.DataLoader
            DataLoader for prediction data
            
        Returns:
        --------
        tuple
            (predictions, nanobody_ids)
        """
        self.model.eval()
        all_preds = []
        all_ids = []
        
        with torch.no_grad():
            for data in loader:
                # Move data to device
                data = data.to(self.device)
                
                # Forward pass
                out = self.model(data)
                
                # Store predictions and IDs
                all_preds.append(out.cpu().numpy())
                
                # Extract nanobody IDs if available
                if hasattr(data, 'nanobody_id'):
                    all_ids.extend(data.nanobody_id)
        
        # Concatenate predictions
        all_preds = np.concatenate(all_preds, axis=0)
        
        return all_preds, all_ids
    
    def predict_3d_structure(self, loader, output_dir='predicted_structures'):
        """
        Predict optimal conditions and generate 3D structures.
        
        Parameters:
        -----------
        loader : torch_geometric.data.DataLoader
            DataLoader for prediction data
        output_dir : str
            Directory to save predicted structures
            
        Returns:
        --------
        dict
            Dictionary of predictions and structure files
        """
        # This method would integrate with a 3D structure prediction tool
        # For now, we'll just make predictions and simulate structure generation
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Make predictions
        predictions, nanobody_ids = self.predict(loader)
        
        results = {}
        
        for i, nanobody_id in enumerate(nanobody_ids):
            pred_ph, pred_ionic = predictions[i]
            
            # Store prediction
            results[nanobody_id] = {
                'predicted_pH': float(pred_ph),
                'predicted_ionic_strength': float(pred_ionic),
                'structure_file': f"{nanobody_id}_predicted.pdb"
            }
            
            # In a real implementation, you would call a structure prediction
            # tool here using the predicted pH and ionic strength
            # For example:
            # predicted_structure = generate_3d_structure(
            #     nanobody_id, pred_ph, pred_ionic
            # )
            
            # Write a placeholder structure file
            with open(os.path.join(output_dir, f"{nanobody_id}_predicted.pdb"), 'w') as f:
                f.write(f"PLACEHOLDER STRUCTURE FOR {nanobody_id}\n")
                f.write(f"PREDICTED pH: {pred_ph:.2f}\n")
                f.write(f"PREDICTED IONIC STRENGTH: {pred_ionic:.2f}\n")
        
        # Save results to JSON
        with open(os.path.join(output_dir, 'predictions.json'), 'w') as f:
            json.dump(results, f, indent=2)
        
        return results
    

class NanobodyFeatureExtractor:
    """
    A class to extract and process nanobody features from graph data.
    """
    
    # Define feature indices with clear names
    FEATURE_INDICES = {
        'atom_type': 0,
        'residue_type': 1,
        'charge': 2,
        'hydrophobicity': 3,
        'h_bond_donor': 4,
        'h_bond_acceptor': 5,
        'aromatic': 6,
        'solvent_accessible': 7,
        'secondary_structure': 8,
        'pka': 9,
        'isoelectric_point': 10,
        # Add more features as needed
    }
    
    @staticmethod
    def get_feature_index(feature_name):
        """Get the index of a feature by name."""
        if feature_name in NanobodyFeatureExtractor.FEATURE_INDICES:
            return NanobodyFeatureExtractor.FEATURE_INDICES[feature_name]
        else:
            raise ValueError(f"Unknown feature name: {feature_name}")
    
    @staticmethod
    def extract_features(x, feature_names):
        """
        Extract features by name from node features.
        
        Parameters:
        -----------
        x : torch.Tensor
            Node features
        feature_names : list
            List of feature names to extract
            
        Returns:
        --------
        dict
            Dictionary of extracted features
        """
        features = {}
        
        for name in feature_names:
            idx = NanobodyFeatureExtractor.get_feature_index(name)
            if idx < x.size(1):
                features[name] = x[:, idx]
            else:
                # Default values if feature doesn't exist
                if name == 'charge':
                    features[name] = torch.zeros(x.size(0))
                elif name == 'hydrophobicity':
                    features[name] = torch.zeros(x.size(0))
                elif name == 'pka':
                    features[name] = torch.ones(x.size(0)) * 7.0
                else:
                    features[name] = torch.zeros(x.size(0))
        
        return features
    
    @staticmethod
    def calculate_aggregate_features(features, batch):
        """
        Calculate aggregate features for each graph in the batch.
        
        Parameters:
        -----------
        features : dict
            Dictionary of node features
        batch : torch.Tensor
            Batch assignment of each node
            
        Returns:
        --------
        dict
            Dictionary of aggregate features
        """
        aggregate_features = {}
        
        # Calculate aggregate features
        for name, feature in features.items():
            # Mean pooling
            agg_mean = global_mean_pool(feature.unsqueeze(1), batch)
            aggregate_features[f'{name}_mean'] = agg_mean.squeeze(-1)
            
            # Sum pooling
            agg_sum = global_add_pool(feature.unsqueeze(1), batch)
            aggregate_features[f'{name}_sum'] = agg_sum.squeeze(-1)
            
            # Min/Max pooling
            if name in ['charge', 'hydrophobicity', 'pka']:
                # Custom aggregation for specific features
                unique_batches = torch.unique(batch)
                agg_min = torch.zeros(len(unique_batches), device=feature.device)
                agg_max = torch.zeros(len(unique_batches), device=feature.device)
                
                for i, b in enumerate(unique_batches):
                    mask = (batch == b)
                    agg_min[i] = torch.min(feature[mask])
                    agg_max[i] = torch.max(feature[mask])
                
                aggregate_features[f'{name}_min'] = agg_min
                aggregate_features[f'{name}_max'] = agg_max
        
        return aggregate_features
    
    @staticmethod
    def calculate_physics_features(features, aggregate_features):
        """
        Calculate physics-based features from extracted features.
        
        Parameters:
        -----------
        features : dict
            Dictionary of node features
        aggregate_features : dict
            Dictionary of aggregate features
            
        Returns:
        --------
        dict
            Dictionary of physics-based features
        """
        physics_features = {}
        
        # Calculate charge-related features
        if 'charge_sum' in aggregate_features:
            # Simplified Debye-Hückel parameter
            physics_features['ionic_strength_factor'] = torch.abs(aggregate_features['charge_sum'])
            
            # Henderson-Hasselbalch factor
            if 'pka_mean' in aggregate_features:
                charge_ratio = torch.sigmoid(aggregate_features['charge_sum']) * 14.0  # Scale to pH range
                physics_features['ph_factor'] = charge_ratio
        
        # Calculate hydrophobicity-related features
        if 'hydrophobicity_sum' in aggregate_features:
            # Simplified hydrophobic effect factor
            physics_features['hydrophobic_factor'] = torch.sigmoid(aggregate_features['hydrophobicity_sum'])
        
        return physics_features
    

class EnhancedNanobodyPINNGNN(nn.Module):
    """
    Enhanced Physics-Informed Graph Neural Network for predicting optimal
    pH and ionic strength for nanobodies.
    """
    
    def __init__(self, node_features, edge_features, hidden_dim=128, num_layers=4, 
                 dropout=0.3, use_gat=True, residual=True, layer_norm=True):
        """
        Initialize the enhanced PINN GNN model.
        
        Parameters:
        -----------
        node_features : int
            Number of node features
        edge_features : int
            Number of edge features
        hidden_dim : int
            Dimension of hidden layers
        num_layers : int
            Number of GNN layers
        dropout : float
            Dropout probability
        use_gat : bool
            Whether to use GAT instead of GCN
        residual : bool
            Whether to use residual connections
        layer_norm : bool
            Whether to use layer normalization
        """
        super(EnhancedNanobodyPINNGNN, self).__init__()
        
        self.node_features = node_features
        self.edge_features = edge_features
        self.hidden_dim = hidden_dim
        self.use_gat = use_gat
        self.residual = residual
        self.layer_norm = layer_norm
        self.feature_extractor = NanobodyFeatureExtractor()
        
        # Node feature transformation
        self.node_encoder = nn.Sequential(
            nn.Linear(node_features, hidden_dim),
            nn.LayerNorm(hidden_dim) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout)
        )
        
        # Edge feature transformation
        self.edge_encoder = nn.Sequential(
            nn.Linear(edge_features, hidden_dim),
            nn.LayerNorm(hidden_dim) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout)
        )
        
        # Graph convolution layers
        self.convs = nn.ModuleList()
        self.layer_norms = nn.ModuleList() if layer_norm else None
        
        # Initialize the first layer with the hidden_dim (after node encoding)
        if use_gat:
            self.convs.append(GATv2Conv(hidden_dim, hidden_dim, edge_dim=hidden_dim))
        else:
            self.convs.append(GCNConv(hidden_dim, hidden_dim))
            
        if layer_norm:
            self.layer_norms.append(nn.LayerNorm(hidden_dim))
        
        # Add additional layers
        for _ in range(num_layers - 1):
            if use_gat:
                self.convs.append(GATv2Conv(hidden_dim, hidden_dim, edge_dim=hidden_dim))
            else:
                self.convs.append(GCNConv(hidden_dim, hidden_dim))
                
            if layer_norm:
                self.layer_norms.append(nn.LayerNorm(hidden_dim))
        
        # Enhanced physics-informed layer
        self.physics_layer = EnhancedPhysicsInformedLayer(hidden_dim, hidden_dim)
        
        # Pooling with attention
        if use_gat:
            # Attention-based pooling
            self.attn_pool = nn.Sequential(
                nn.Linear(hidden_dim, 1),
                nn.Sigmoid()
            )
        else:
            self.attn_pool = None
            
        # Default pooling as fallback
        self.pool = global_mean_pool
        
        # Fully connected layers for pH prediction with residual connections
        self.ph_branch = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        # Fully connected layers for ionic strength prediction with residual connections
        self.ionic_branch = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        # Physics-based constraint modules
        # Henderson-Hasselbalch inspired module for pH
        self.ph_physics = nn.Sequential(
            nn.Linear(hidden_dim + 3, hidden_dim // 2),  # +3 for extra physics features
            nn.LayerNorm(hidden_dim // 2) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        # Debye-Hückel inspired module for ionic strength
        self.ionic_physics = nn.Sequential(
            nn.Linear(hidden_dim + 3, hidden_dim // 2),  # +3 for extra physics features
            nn.LayerNorm(hidden_dim // 2) if layer_norm else nn.Identity(),
            nn.SiLU(),
            nn.Linear(hidden_dim // 2, 1)
        )
        
        # Blending parameters with proper initialization
        self.alpha_ph = nn.Parameter(torch.tensor([0.7]))  # Start with higher weight on NN
        self.alpha_ionic = nn.Parameter(torch.tensor([0.7]))
        
        # Physics contribution scaling (learnable)
        self.physics_scale_ph = nn.Parameter(torch.tensor([0.2]))
        self.physics_scale_ionic = nn.Parameter(torch.tensor([0.3]))
        
        # Batch normalization for final outputs
        self.bn_ph = nn.BatchNorm1d(1)
        self.bn_ionic = nn.BatchNorm1d(1)
        
        # Dropout
        self.dropout = nn.Dropout(dropout)
    
    def attention_pooling(self, x, batch):
        """Custom attention-based pooling."""
        if self.attn_pool is None:
            return self.pool(x, batch)
            
        # Calculate attention weights
        attention = self.attn_pool(x)
        
        # Apply attention weights
        x_weighted = x * attention
        
        # Pool with attention weights
        return global_add_pool(x_weighted, batch)
    
    def forward(self, data):
        """
        Forward pass through the model.
        
        Parameters:
        -----------
        data : torch_geometric.data.Data
            Input graph data
            
        Returns:
        --------
        torch.Tensor
            Predicted [pH, ionic_strength]
        """
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        
        # Extract features using the feature extractor
        features = self.feature_extractor.extract_features(
            x, ['charge', 'hydrophobicity', 'pka']
        )
        
        # Transform node and edge features
        x = self.node_encoder(x)
        edge_attr = self.edge_encoder(edge_attr)
        
        # Store original x for residual connections
        x_prev = x
        
        # Apply GNN layers
        for i, conv in enumerate(self.convs):
            if self.use_gat:
                x_new = conv(x, edge_index, edge_attr)
            else:
                x_new = conv(x, edge_index)
            
            # Apply layer normalization if enabled
            if self.layer_norm:
                x_new = self.layer_norms[i](x_new)
            
            # Apply activation and dropout
            x_new = F.silu(x_new)
            x_new = self.dropout(x_new)
            
            # Add residual connection if enabled
            if self.residual and i > 0:  # Skip first layer for dim matching
                x = x_new + x
            else:
                x = x_new
        
        # Apply enhanced physics-informed layer
        x = self.physics_layer(
            x, 
            features['charge'], 
            features['hydrophobicity'],
            features.get('pka', None)
        )
        
        # Calculate aggregate features
        agg_features = self.feature_extractor.calculate_aggregate_features(features, batch)
        
        # Calculate physics-based features
        physics_features = self.feature_extractor.calculate_physics_features(features, agg_features)
        
        # Global pooling with attention
        x_pooled = self.attention_pooling(x, batch)
        
        # Prepare physics inputs - concatenate with relevant physics features
        physics_inputs = []
        
        # Add relevant aggregate features
        for key in ['charge_sum', 'hydrophobicity_mean', 'pka_mean']:
            if key in agg_features:
                physics_inputs.append(agg_features[key].unsqueeze(1))
            else:
                # Use zeros as fallback
                physics_inputs.append(torch.zeros_like(x_pooled[:, :1]))
                
        # Concatenate pooled features with physics inputs
        physics_concat = torch.cat([x_pooled] + physics_inputs, dim=1)
        
        # Neural network pathways
        ph_pred_nn = self.ph_branch(x_pooled)
        ionic_pred_nn = self.ionic_branch(x_pooled)
        
        # Physics-informed pathways
        ph_pred_phys = self.ph_physics(physics_concat)
        ionic_pred_phys = self.ionic_physics(physics_concat)
        
        # Apply physics constraints
        ph_constraint = torch.zeros_like(ph_pred_nn)
        ionic_constraint = torch.zeros_like(ionic_pred_nn)
        
        if 'ph_factor' in physics_features:
            ph_constraint = physics_features['ph_factor'].unsqueeze(1)
            
        if 'ionic_strength_factor' in physics_features:
            ionic_constraint = physics_features['ionic_strength_factor'].unsqueeze(1)
        
        # Blend neural network and physics-based predictions
        alpha_ph = torch.sigmoid(self.alpha_ph)  # Ensure between 0 and 1
        alpha_ionic = torch.sigmoid(self.alpha_ionic)
        
        # Calculate physics scaling factors
        physics_scale_ph = F.softplus(self.physics_scale_ph)
        physics_scale_ionic = F.softplus(self.physics_scale_ionic)
        
        # Blend predictions
        ph_pred = alpha_ph * ph_pred_nn + (1 - alpha_ph) * ph_pred_phys + physics_scale_ph * ph_constraint
        ionic_pred = alpha_ionic * ionic_pred_nn + (1 - alpha_ionic) * ionic_pred_phys + physics_scale_ionic * ionic_constraint
        
        # Apply batch normalization
        ph_pred = self.bn_ph(ph_pred)
        ionic_pred = self.bn_ionic(ionic_pred)
        
        # Ensure predictions are in valid ranges
        # pH typically ranges from 0 to 14
        ph_pred = 7.0 + torch.tanh(ph_pred) * 7.0  # Center around 7 with range 0-14
        
        # Ionic strength is non-negative (typically 0 to 1 mol/L for most buffers)
        ionic_pred = F.softplus(ionic_pred)  # Ensures positive values
        
        # Combine predictions
        predictions = torch.cat([ph_pred, ionic_pred], dim=1)
        
        return predictions