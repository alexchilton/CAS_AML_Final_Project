import os
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
from torch_geometric.data import Data, Dataset
import torch.nn.functional as F


def process_protein_with_chains(pdb_id, pdb_dir):
    """Process a protein with multiple chains."""
    # Find all chain-specific data files
    chain_files = {}
    for filename in os.listdir(pdb_dir):
        if filename.endswith("_data.pkl") and not filename.startswith(pdb_id + "_"):
            # This is likely a chain-specific file
            chain_id = filename.replace("_data.pkl", "").split("_")[-1]
            chain_files[chain_id] = os.path.join(pdb_dir, filename)
    
    if not chain_files:
        # If no chain-specific files, use the main data file
        data_path = os.path.join(pdb_dir, f"{pdb_id}_data.pkl")
        if os.path.exists(data_path):
            with open(data_path, 'rb') as f:
                pdb_data = pickle.load(f)
            
            # Extract chains from the data
            chains = pdb_data.get('residues_by_chain', {})
            return {chain_id: pdb_data for chain_id in chains}
    
    # Load data for each chain
    chain_data = {}
    for chain_id, filepath in chain_files.items():
        with open(filepath, 'rb') as f:
            chain_data[chain_id] = pickle.load(f)
    
    return chain_data

def extract_chemical_features(pdb_data, chain_id):
    """Extract relevant chemical features that could influence pH optimum."""
    features = []
    
    # 1. Number of charged residues (Asp, Glu, Lys, Arg, His)
    residues = pdb_data.get('residues_by_chain', {}).get(chain_id, [])
    acid_count = sum(1 for _, res_name in residues if res_name in ['ASP', 'GLU'])
    base_count = sum(1 for _, res_name in residues if res_name in ['LYS', 'ARG'])
    his_count = sum(1 for _, res_name in residues if res_name == 'HIS')
    
    # 2. Surface exposure of charged residues (approximation)
    # We could improve this by calculating solvent accessible surface area
    
    # 3. Charge distribution
    charges = pdb_data.get('charges', {}).get(chain_id, {})
    total_charge = sum(charges.values())
    charge_density = total_charge / len(residues) if residues else 0
    
    # 4. Secondary structure composition
    ss_data = pdb_data.get('secondary_structure', {}).get(chain_id, {})
    helix_count = sum(1 for ss in ss_data.values() if ss == 'H')
    sheet_count = sum(1 for ss in ss_data.values() if ss == 'E')
    helix_percent = helix_count / len(ss_data) if ss_data else 0
    sheet_percent = sheet_count / len(ss_data) if ss_data else 0
    
    # 5. Hydrophobic residues
    hydrophobic = pdb_data.get('hydrophobic_residues', {}).get(chain_id, set())
    hydrophobic_percent = len(hydrophobic) / len(residues) if residues else 0
    
    # 6. Torsion angle distributions (might indicate structural flexibility)
    torsions = pdb_data.get('torsions', {}).get(chain_id, {})
    torsion_mean = np.mean(list(torsions.values())) if torsions else 0
    torsion_std = np.std(list(torsions.values())) if torsions else 0
    
    # Combine features
    features = [
        acid_count, base_count, his_count,
        acid_count / (len(residues) if residues else 1),  # Acid ratio
        base_count / (len(residues) if residues else 1),  # Base ratio
        his_count / (len(residues) if residues else 1),   # His ratio
        total_charge, charge_density,
        helix_percent, sheet_percent,
        hydrophobic_percent,
        torsion_mean, torsion_std
    ]
    
    return features

def analyze_feature_importance(model, loader, device):
    """Analyze feature importance using permutation importance."""
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.inspection import permutation_importance
    
    # Get predictions and features
    model.eval()
    all_preds = []
    all_features = []
    all_targets = []
    
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            
            # Get node-level features
            node_features = data.x.cpu().numpy()
            
            # Forward pass
            pred = model(data.x, data.edge_index, 
                         data.edge_attr if hasattr(data, 'edge_attr') else None, 
                         data.batch)
            
            # Store prediction and target
            all_preds.extend(pred.cpu().numpy())
            all_targets.extend(data.y.cpu().numpy())
            
            # For each graph in the batch, get mean feature values
            batch_indices = data.batch.cpu().numpy()
            for i in range(data.num_graphs):
                # Get nodes for this graph
                mask = batch_indices == i
                graph_features = node_features[mask]
                
                # Calculate mean feature values for this graph
                mean_features = np.mean(graph_features, axis=0)
                all_features.append(mean_features)
    
    # Convert to numpy arrays
    X = np.array(all_features)
    y = np.array(all_targets)
    
    # Define feature names - adjust based on your actual feature dimensions
    feature_names = [
        "AA_ALA", "AA_ARG", "AA_ASN", "AA_ASP", "AA_CYS", 
        "AA_GLN", "AA_GLU", "AA_GLY", "AA_HIS", "AA_ILE", 
        "AA_LEU", "AA_LYS", "AA_MET", "AA_PHE", "AA_PRO", 
        "AA_SER", "AA_THR", "AA_TRP", "AA_TYR", "AA_VAL",
        "SS_H", "SS_E", "SS_C", "SS_B", "SS_G", "SS_I", "SS_T", "SS_S",
        "Hydrophobicity", "Charge"
    ]
    
    # Calculate feature importance using permutation importance
    # First train a simple model on the extracted features
    from sklearn.ensemble import RandomForestRegressor
    rf = RandomForestRegressor(n_estimators=100, random_state=42)
    rf.fit(X, y)
    
    # Calculate permutation importance
    result = permutation_importance(
        rf, X, y, n_repeats=10, random_state=42, n_jobs=-1
    )
    
    # Get importance scores
    importance_scores = result.importances_mean
    
    # Print correlations between features and predictions
    for i, feature_name in enumerate(feature_names[:X.shape[1]]):  # Ensure we only use available features
        if i < X.shape[1]:  # Check if the feature index is valid
            feature_values = X[:, i]
            correlation = np.corrcoef(feature_values, y.flatten())[0, 1]
            print(f"Correlation of {feature_name} with pH: {correlation:.4f}")
    
    # Sort features by importance
    indices = np.argsort(importance_scores)[::-1]
    sorted_feature_names = [feature_names[i] for i in indices if i < X.shape[1]]
    sorted_importance = importance_scores[indices]
    
    # Create importance plot
    plt.figure(figsize=(12, 8))
    plt.bar(range(len(sorted_importance)), sorted_importance)
    plt.xticks(range(len(sorted_importance)), sorted_feature_names, rotation=90)
    plt.xlabel('Features')
    plt.ylabel('Importance Score (Mean decrease in impurity)')
    plt.title('Feature Importance for pH Prediction')
    plt.tight_layout()
    plt.savefig('feature_importance.png')
    plt.show()
    
    return importance_scores, feature_names[:X.shape[1]]