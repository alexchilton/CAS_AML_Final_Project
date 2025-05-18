import torch
from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from dataloader import  ProteinGraphDataset
from model import GNN

# Prepare dataset
processed_dir_path='C:/Users/laran/Documents/code/dataset/prot_selection/enzyme_ph_structures/new_protein_data_networkx'
dataset = ProteinGraphDataset(root='data', processed_dir=processed_dir_path)

# Split dataset into train/val/test
train_indices, test_indices = train_test_split(
    range(len(dataset)), test_size=0.2, random_state=42)
train_indices, val_indices = train_test_split(
    train_indices, test_size=0.25, random_state=42)  # 0.25 * 0.8 = 0.2 of total

train_dataset = [dataset[i] for i in train_indices]
val_dataset = [dataset[i] for i in val_indices]
test_dataset = [dataset[i] for i in test_indices]

# Create data loaders
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# Calculate number of features
sample_data = dataset[0]
num_node_features = sample_data.x.size(1)
num_edge_features = sample_data.edge_attr.size(1) if hasattr(sample_data, 'edge_attr') else 0

# Initialize model
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = GNN(num_node_features, num_edge_features).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.001, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5)
criterion = torch.nn.MSELoss()

# Training function
def train():
    model.train()
    total_loss = 0
    
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        
        # Forward pass
        out = model(data.x, data.edge_index, 
                    data.edge_attr if hasattr(data, 'edge_attr') else None, 
                    data.batch)
        
        # Compute loss
        loss = criterion(out.squeeze(), data.y)
        
        # Backward pass
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item() * data.num_graphs
    
    return total_loss / len(train_dataset)

# Evaluation function
def evaluate(loader):
    model.eval()
    total_loss = 0
    predictions = []
    actual = []
    
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            
            # Forward pass
            out = model(data.x, data.edge_index, 
                        data.edge_attr if hasattr(data, 'edge_attr') else None, 
                        data.batch)
            
            # Compute loss
            loss = criterion(out.squeeze(), data.y)
            total_loss += loss.item() * data.num_graphs
            
            # Store predictions and actual values
            predictions.extend(out.squeeze().cpu().numpy())
            actual.extend(data.y.cpu().numpy())
    
    # Calculate metrics
    mse = mean_squared_error(actual, predictions)
    mae = mean_absolute_error(actual, predictions)
    r2 = r2_score(actual, predictions)
    
    return total_loss / len(loader.dataset), mse, mae, r2, predictions, actual

# Training loop
best_val_loss = float('inf')
patience = 15
patience_counter = 0
train_losses = []
val_losses = []

for epoch in range(1, 201):
    # Train
    train_loss = train()
    train_losses.append(train_loss)
    
    # Validate
    val_loss, val_mse, val_mae, val_r2, _, _ = evaluate(val_loader)
    val_losses.append(val_loss)
    
    # Print progress
    print(f'Epoch: {epoch:03d}, Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}, '
          f'Val MSE: {val_mse:.4f}, Val MAE: {val_mae:.4f}, Val R2: {val_r2:.4f}')
    
    # Update learning rate
    scheduler.step(val_loss)
    
    # Early stopping
    if val_loss < best_val_loss:
        best_val_loss = val_loss
        patience_counter = 0
        # Save best model
        torch.save(model.state_dict(), 'best_ph_model.pt')
    else:
        patience_counter += 1
        if patience_counter >= patience:
            print(f'Early stopping at epoch {epoch}')
            break

# Plot training curves
plt.figure(figsize=(10, 6))
plt.plot(train_losses, label='Train Loss')
plt.plot(val_losses, label='Validation Loss')
plt.xlabel('Epoch')
plt.ylabel('Loss (MSE)')
plt.title('Training and Validation Loss')
plt.legend()
plt.savefig('training_curve.png')
plt.show()

# Load best model and evaluate on test set
model.load_state_dict(torch.load('best_ph_model.pt'))
test_loss, test_mse, test_mae, test_r2, test_preds, test_actual = evaluate(test_loader)

print(f'Test Results - Loss: {test_loss:.4f}, MSE: {test_mse:.4f}, '
      f'MAE: {test_mae:.4f}, R2: {test_r2:.4f}')

# Plot predictions vs actual values
plt.figure(figsize=(10, 6))
plt.scatter(test_actual, test_preds, alpha=0.6)
plt.plot([min(test_actual), max(test_actual)], [min(test_actual), max(test_actual)], 'r--')
plt.xlabel('Actual pH Optimum')
plt.ylabel('Predicted pH Optimum')
plt.title('Predicted vs Actual pH Optimum Values')
plt.savefig('prediction_scatter.png')
plt.show()