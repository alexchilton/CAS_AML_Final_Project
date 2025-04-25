# test script for self-supervised learning
import os
import torch
import numpy as np
import torch.nn.functional as F
from torch_geometric.data import DataLoader

# Import your existing modules
from nanobody_graph_builder import NanobodyGraphBuilder
from nanobody_pinn_gnn import NanobodyPINNGNN, NanobodyDataset, NanobodyPINNTrainer

# Create a custom self-supervised trainer that doesn't require labels
class SelfSupervisedTrainer:
    """Self-supervised trainer for nanobody PINN GNN model."""
    
    def __init__(self, model, optimizer, device=None):
        self.model = model
        self.optimizer = optimizer
        self.device = device if device is not None else torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        self.history = {'train_loss': [], 'val_loss': []}
    
    def self_supervised_loss(self, out, data):
        """Calculate loss without ground truth labels."""
        # Extract pH and ionic strength predictions
        ph_pred = out[:, 0]
        ionic_pred = out[:, 1]
        
        # Extract physics-relevant features from the data if available
        # (You'll need to adapt this based on how your features are stored)
        if hasattr(data, 'x'):
            # Extract charge information (assuming it's in a specific position)
            # Adjust these indices based on your actual feature ordering
            charge_idx = -2  # Example: charge is second-to-last feature
            hydrophobic_idx = -3  # Example: hydrophobicity is third-to-last feature
            
            # Calculate aggregate charge and hydrophobicity
            charge = data.x[:, charge_idx]
            hydrophobic = data.x[:, hydrophobic_idx]
            
            # Aggregate by taking mean across nodes
            avg_charge = torch.mean(charge)
            avg_hydrophobic = torch.mean(hydrophobic)
            
            # Physics-based constraints:
            
            # 1. Henderson-Hasselbalch inspired constraint
            # Higher negative charge generally correlates with higher pH
            # pH ~ 7 + scaled_charge
            scaled_charge = torch.sigmoid(avg_charge) * 14.0 - 7.0  # Scale to pH range
            expected_ph = 7.0 + scaled_charge
            ph_loss = F.mse_loss(ph_pred, expected_ph.expand_as(ph_pred))
            
            # 2. Ionic strength constraint
            # Higher absolute charge generally requires higher ionic strength
            abs_charge = torch.abs(avg_charge)
            expected_ionic = 0.1 + 0.4 * torch.sigmoid(abs_charge)  # Scale to reasonable range
            ionic_loss = F.mse_loss(ionic_pred, expected_ionic.expand_as(ionic_pred))
            
            # 3. pH stability constraint (most proteins stable around pH 7)
            stability_loss = torch.mean(torch.abs(ph_pred - 7.0))
            
            # Combine losses
            total_loss = ph_loss + ionic_loss + 0.5 * stability_loss
        else:
            # Fallback if no features are available - just use reasonable range constraints
            ph_range_loss = torch.mean(F.relu(torch.abs(ph_pred - 7.0) - 2.0))  # Penalize pH outside 5-9
            ionic_range_loss = torch.mean(F.relu(ionic_pred - 0.5) + F.relu(0.05 - ionic_pred))  # Penalize outside 0.05-0.5
            total_loss = ph_range_loss + ionic_range_loss
        
        return total_loss
    
    def train_epoch(self, train_loader):
        self.model.train()
        total_loss = 0
        
        for data in train_loader:
            data = data.to(self.device)
            self.optimizer.zero_grad()
            out = self.model(data)
            loss = self.self_supervised_loss(out, data)
            loss.backward()
            self.optimizer.step()
            total_loss += loss.item() * data.num_graphs
            
        return total_loss / len(train_loader.dataset)
    
    def evaluate(self, loader):
        self.model.eval()
        total_loss = 0
        
        with torch.no_grad():
            for data in loader:
                data = data.to(self.device)
                out = self.model(data)
                loss = self.self_supervised_loss(out, data)
                total_loss += loss.item() * data.num_graphs
                
        return total_loss / len(loader.dataset)
    
    def train(self, train_loader, val_loader, num_epochs=10, checkpoint_dir='checkpoints'):
        os.makedirs(checkpoint_dir, exist_ok=True)
        best_val_loss = float('inf')
        
        for epoch in range(1, num_epochs + 1):
            train_loss = self.train_epoch(train_loader)
            val_loss = self.evaluate(val_loader)
            
            self.history['train_loss'].append(train_loss)
            self.history['val_loss'].append(val_loss)
            
            print(f'Epoch {epoch}/{num_epochs}, Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}')
            
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save({
                    'epoch': epoch,
                    'model_state_dict': self.model.state_dict(),
                    'optimizer_state_dict': self.optimizer.state_dict(),
                    'val_loss': val_loss
                }, os.path.join(checkpoint_dir, 'best_model.pth'))
                print(f'Model saved with validation loss: {val_loss:.4f}')
                
        return self.history

# Main test function
def test_self_supervised_model():
    # Parameters
    data_dir = "nanos_networkx"
    batch_size = 4
    hidden_dim = 64
    num_epochs = 5
    
    # Define the pytorch_geometric directory
    dataset_dir = os.path.join(data_dir, "pytorch_geometric")
    
    # Check if dataset exists
    if os.path.exists(dataset_dir) and os.path.exists(os.path.join(dataset_dir, 'splits.json')):
        print(f"Using existing dataset at {dataset_dir}")
        # Load splits from file
        import json
        with open(os.path.join(dataset_dir, 'splits.json'), 'r') as f:
            splits_info = json.load(f)
            
        # Take just a few samples from each split
        max_samples = 10
        train_ids = splits_info['train'][:min(max_samples, len(splits_info['train']))]
        val_ids = splits_info['val'][:min(max_samples//2, len(splits_info['val']))]
        test_ids = splits_info['test'][:min(max_samples//2, len(splits_info['test']))]
    else:
        print("Dataset not found. Please create it first.")
        return None, None
    
    # Create file paths
    train_paths = [os.path.join(dataset_dir, f"{nanobody_id}.pt") for nanobody_id in train_ids]
    val_paths = [os.path.join(dataset_dir, f"{nanobody_id}.pt") for nanobody_id in val_ids]
    test_paths = [os.path.join(dataset_dir, f"{nanobody_id}.pt") for nanobody_id in test_ids]
    
    # Create datasets
    train_dataset = NanobodyDataset(train_paths)
    val_dataset = NanobodyDataset(val_paths)
    test_dataset = NanobodyDataset(test_paths)
    
    print(f"Dataset sizes: Train {len(train_dataset)}, Val {len(val_dataset)}, Test {len(test_dataset)}")
    
    # Create data loaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    # Load a sample to get feature dimensions
    print("Loading a sample to get feature dimensions...")
    sample_data = train_dataset[0]
    node_features = sample_data.x.size(1)
    edge_features = sample_data.edge_attr.size(1)
    
    print(f"Sample nanobody has {sample_data.num_nodes} nodes, {sample_data.num_edges} edges")
    print(f"Node feature dimension: {node_features}, Edge feature dimension: {edge_features}")
    
    # Create model
    model = NanobodyPINNGNN(
        node_features=node_features,
        edge_features=edge_features,
        hidden_dim=hidden_dim,
        num_layers=3,
        dropout=0.2,
        use_gat=True
    )
    
    # Create optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    
    # Create self-supervised trainer
    trainer = SelfSupervisedTrainer(
        model=model,
        optimizer=optimizer
    )
    
    # Create output directory for checkpoints
    checkpoint_dir = "test_checkpoints"
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    # Train for a few epochs
    print("Starting self-supervised training for a few epochs...")
    history = trainer.train(
        train_loader=train_loader,
        val_loader=val_loader,
        num_epochs=num_epochs
    )
    
    # Print final metrics
    print("\nTraining complete!")
    print(f"Final train loss: {history['train_loss'][-1]:.4f}")
    print(f"Final validation loss: {history['val_loss'][-1]:.4f}")
    
    return trainer, history

if __name__ == "__main__":
    print("Testing self-supervised Nanobody PINN GNN model")
    trainer, history = test_self_supervised_model()
    print("Test completed successfully")