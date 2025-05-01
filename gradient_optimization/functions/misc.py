from sklearn.manifold import MDS
import numpy as np
import torch
import copy

import matplotlib.pyplot as plt

def compare_images(cmaps: list, titles: list, masking=False, masks=None):
      
    fig, axes = plt.subplots(1, len(cmaps), figsize=(len(cmaps)*3.5, 3.5))  

    assert len(titles)==len(cmaps), 'Lists not same length'

    if masking:
        cmaps = copy.deepcopy(cmaps)
        for i, contact_map in enumerate(cmaps):
            contact_map[~masks[i]] = 0.5

    for i in range(len(cmaps)):
        axes[i].imshow(cmaps[i], cmap='viridis')
        axes[i].set_title(titles[i])
        axes[i].axis('off')
        
    plt.tight_layout()
    plt.show()

def get_range_matrix(contact_maps:list, mask, configs):

    stacked = np.stack(contact_maps, axis=2)
    range_matrix = np.zeros((contact_maps[0].shape[0], contact_maps[0].shape[1], 2))
    
    for i in range(contact_maps[0].shape[0]):
        for ii in range(contact_maps[0].shape[1]):
            if np.sum(stacked[i, ii])==0: continue
            idx = stacked[i, ii].argmax()
            range_matrix[i,ii] = np.array([configs[idx]['lower'], configs[idx]['upper']])
            range_matrix[ii,i] = np.array([configs[idx]['lower'], configs[idx]['upper']])

    return range_matrix

def get_cmaps_from_range_matrix(range_matrix, configs):
    N = range_matrix.shape[0]
    num_configs = len(configs)
    contact_maps = [np.zeros((N, N), dtype=np.float32) for _ in range(num_configs)]

    for i in range(N):
        for j in range(N):
            lower, upper = range_matrix[i, j]
            if lower == 0 and upper == 0:
                continue  # skip unknowns or empty

            for idx, config in enumerate(configs):
                if np.isclose(lower, config['lower']) and np.isclose(upper, config['upper']):
                    contact_maps[idx][i, j] = 1.0
                    break  # stop after first match (only one config per pair)

    return contact_maps

def fill_missing(distance_ranges, known_mask, configs, max_passes=10, verbose=True):
    N = distance_ranges.shape[0]

    for p in range(max_passes):
        updates = 0

        for i in range(N):
            for j in range(i+1, N):  # Only upper triangle
                if known_mask[i, j]:
                    continue  # Already filled

                candidate_ranges = []
                for k in range(N):
                    if known_mask[i, k] and known_mask[k, j]:
                        min_ik, max_ik = distance_ranges[i, k]
                        min_kj, max_kj = distance_ranges[k, j]

                        lower_sum = min_ik + min_kj
                        upper_sum = max_ik + max_kj

                        if lower_sum <= upper_sum:
                            candidate_ranges.append((lower_sum, upper_sum))

                if candidate_ranges:
                    # Determine overall lower and upper bounds
                    inferred_lower = min([c[0] for c in candidate_ranges])
                    inferred_upper = max([c[1] for c in candidate_ranges])

                    # Find all configs that contain [inferred_lower, inferred_upper]
                    valid_configs = []
                    for config in configs:
                        if config['lower'] <= inferred_lower and config['upper'] >= inferred_upper:
                            valid_configs.append(config)

                    if valid_configs:
                        # Pick config with minimal (upper - lower) width
                        selected_config = min(valid_configs, key=lambda x: (x['upper'] - x['lower']))

                        # Fill
                        distance_ranges[i, j] = np.array([selected_config['lower'], selected_config['upper']])
                        distance_ranges[j, i] = np.array([selected_config['lower'], selected_config['upper']])
                        known_mask[i, j] = True
                        known_mask[j, i] = True
                        updates += 1

        if verbose:
            print(f"Pass {p}: {updates} updates")
        if updates == 0:
            break

    return distance_ranges

def reconstruct_coords_local_plural_maps(contact_maps: list, configs: list, hard_mask=None, sharpness=10., window=25, dim=3, lr=1e-3, max_iter=1000, prints=10):
    """
    Reconstruct 3D coordinates from multiple contact maps with different distance constraints.
    
    Args:
        contact_maps: List of binary contact maps
        configs: List of dictionaries containing 'lower' and 'upper' bounds for each map
        window: Window size for local interactions
        dim: Dimensionality of the coordinate space (default: 3D)
        lr: Learning rate for optimization
        max_iter: Maximum number of iterations
        prints: Number of progress updates
        
    Returns:
        Reconstructed coordinates as numpy array
    """

    print ('Gradient-based optimization with multiple contact maps...')
    
    min_distance = 2.3  # Minimum allowed distance between points
    
    N = contact_maps[0].shape[0]
    # Ensure all contact maps have the same shape
    assert len(set([cmap.shape for cmap in contact_maps])) == 1, "Contact maps must have the same dimensions"
    
    # Initialize coordinates as a linear chain to provide better starting point
    coords = advanced_initialize_coords(N, dim, contact_maps, configs)
    coords.requires_grad = True
    
    # Use a more robust optimizer and a learning rate scheduler
    optimizer = torch.optim.Adam([coords], lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=50)
    
    # Convert contact maps to tensors
    tensor_maps = [torch.tensor(contact_map, dtype=torch.float32) for contact_map in contact_maps]
    
    # Create repulsion mask (all non-diagonal elements)
    repulsion_mask = (torch.eye(N) == 0)

    if hard_mask is None:
        # Create a mask for local interactions
        local_mask = torch.zeros_like(tensor_maps[0])
        for i in range(N):
            for j in range(max(0, i - window), min(N, i + window + 1)):
                local_mask[i, j] = 1
    else:
        local_mask = torch.tensor(hard_mask.astype(bool))
    
    # Calculate positive weights for each contact map - ensure they're positive
    pos_weights = []
    for i, _ in enumerate(configs):
        # Make sure pos_weight is positive by using absolute value or a minimum bound
        ratio = max(1.0, local_mask.sum() / (tensor_maps[i].sum() + 1e-8))
        pos_weight = (ratio - 1.0).detach()
        pos_weight = max(1.0, pos_weight)  # Ensure positive weight
        pos_weights.append(pos_weight)
    
    # Best loss tracking for early stopping
    best_loss = float('inf')
    patience_counter = 0
    best_coords = coords.clone().detach()
    
    # Improved BCE loss function
    bce_loss_fn = torch.nn.BCEWithLogitsLoss(reduction='none')
    
    for step in range(max_iter):
        optimizer.zero_grad()
        
        # Calculate pairwise distances
        dists = torch.cdist(coords, coords, p=2)
        
        # Sequential connectivity loss - ensuring connected residues are close
        sequential_dists = torch.diag(dists, diagonal=1)
        connectivity_loss = torch.mean((sequential_dists - 3.8) ** 2)  # ~3.8Å is typical C-alpha distance
        
        # Calculate contact predictions for each config
        total_loss = torch.tensor(0.0, requires_grad=True)
        partial_losses = []
        
        for i, config in enumerate(configs):
            # Calculate logits for BCE loss
            if config['lower'] == 0:
                # For upper bound only, predict contact if distance < upper
                logits = sharpness * (config['upper'] - dists) * 0.1 * torch.rand(dists.shape)
            else:
                # For range prediction, use a composite sigmoid
                upper_logits = sharpness * (config['upper'] - dists)
                lower_logits = sharpness * (dists - config['lower'])
                # Combine logits - both conditions must be true
                logits = torch.min(upper_logits, lower_logits)
            
            # Apply mask to focus on meaningful interactions
            weight_matrix = torch.ones_like(tensor_maps[i])
            weight_matrix = weight_matrix * pos_weights[i] * tensor_maps[i] + (1 - tensor_maps[i])
            
            # Calculate loss with weights
            map_loss = bce_loss_fn(logits, tensor_maps[i])
            weighted_loss = (map_loss * weight_matrix).mean() # original .mean()
            partial_losses.append(weighted_loss)
            total_loss = total_loss + weighted_loss
            
        # Add repulsion loss to maintain minimum distances
        min_dist_violation = torch.relu(min_distance - dists)
        repulsion_loss = (min_dist_violation[repulsion_mask] ** 2).sum()
        
        # Add all regularization terms
        total_loss = total_loss + 0.1 * repulsion_loss + 0.05 * connectivity_loss
        
        # Backwards pass and optimization
        total_loss.backward()
        
        # Gradient clipping to prevent instability
        torch.nn.utils.clip_grad_norm_([coords], max_norm=1.0)
        
        optimizer.step()
        scheduler.step(total_loss)
        
        # Print progress and update best model if needed
        if step % max(1, int(max_iter // prints)) == 0:
            partial_losses_str = [f"{p_loss.item():8.2f}" for p_loss in partial_losses]
            print(f"Step {step:7d} | total loss: {total_loss.item():8.2f} ({partial_losses_str})")
            
            if total_loss.item() < best_loss:
                best_loss = total_loss.item()
                best_coords = coords.clone().detach()
                patience_counter = 0
            else:
                patience_counter += 1
                
            # Early stopping
            if patience_counter > 100:
                print("Early stopping triggered.")
                break
    
    # Final loss report
    partial_losses_str = [f"{p_loss.item():8.2f}" for p_loss in partial_losses]
    print(f"Step {step:7d} | total loss: {total_loss.item():8.2f} ({partial_losses_str}), final.")
    print ('')
    
    # Return best coordinates
    return best_coords.detach().numpy()

def advanced_initialize_coords(N, dim=3, contact_maps=None, configs=None):
    coords = torch.zeros((N, dim))
    
    # Si des cartes de contact sont disponibles, utilisez-les pour l'initialisation
    if contact_maps is not None:
        
        # Créez une matrice de distance approximative à partir des contacts
        dist_mat = np.ones((N, N)) * 15.0  # Distance par défaut grande
        for i in range(N):
            for j in range(N):
                if contact_maps[0][i, j] > 0.5:  # S'il y a un contact
                    dist_mat[i, j] = configs[0]['upper'] * 0.8  # Un peu moins que la limite supérieure
        
        # Appliquez MDS pour obtenir les coordonnées initiales
        mds = MDS(n_components=dim, dissimilarity='precomputed', random_state=42)
        coords_np = mds.fit_transform(dist_mat)
        coords = torch.tensor(coords_np, dtype=torch.float32)
    else:
        # Sinon, utilisez l'initialisation zigzag améliorée
        for i in range(N):
            coords[i, 0] = i * 3.8
            if dim > 1:
                coords[i, 1] = 2.0 * np.sin(i * 0.5)
            if dim > 2:
                coords[i, 2] = 2.0 * np.cos(i * 0.5)
    
    return coords


def reconstruct_coords_from_contact_map(contact_map, contact_dist=8.0, non_contact_dist=20.0):
    N = contact_map.shape[0]
    
    dist_matrix = np.where(contact_map == 1, contact_dist, non_contact_dist)
    np.fill_diagonal(dist_matrix, 0.0)  # Distance to self is 0

    mds = MDS(n_components=3, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(dist_matrix)

    return coords
