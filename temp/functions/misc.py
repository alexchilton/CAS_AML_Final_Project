from sklearn.manifold import MDS
import numpy as np
import torch
import copy

from collections import defaultdict
import heapq
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt

def compute_dRMSD(coords1, coords2):
    dists1 = squareform(pdist(coords1))
    dists2 = squareform(pdist(coords2))
    
    diff = dists1 - dists2
    N = coords1.shape[0]
    return np.sqrt(np.sum(diff**2) / (N*(N-1)))

def kabsch(P, Q):
    """Optimal rotation using Kabsch algorithm."""
    C = np.dot(np.transpose(P), Q)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    return U

def superimpose(P, Q):
    """Align P onto Q using Kabsch and return aligned P."""
    P_centroid = P.mean(axis=0)
    Q_centroid = Q.mean(axis=0)
    P_centered = P - P_centroid
    Q_centered = Q - Q_centroid
    U = kabsch(P_centered, Q_centered)
    return np.dot(P_centered, U) + Q_centroid

def compute_rmsd(P, Q):
    """Compute RMSD after optimal alignment."""
    P_aligned = superimpose(P, Q)
    score = np.sqrt(np.mean(np.sum((P_aligned - Q)**2, axis=1)))
    return score

def compute_tmscore(P, Q):
    """Simple TM-score approximation from aligned coordinates."""
    d0 = 1.24 * (len(P) - 15) ** (1 / 3) - 1.8
    P_aligned = superimpose(P, Q)
    dist2 = np.sum((P_aligned - Q)**2, axis=1)
    score = np.sum(1 / (1 + dist2 / (d0**2))) / len(P)
    return score

def compare_images(cmaps: list, titles: list, masking=False, masks=None, suptitle=None, save=False):
      
    fig, axes = plt.subplots(1, len(cmaps), figsize=(len(cmaps)*3.5, 3.5))  

    assert len(titles)==len(cmaps), 'Lists not same length'

    if masking:
        cmaps = copy.deepcopy(cmaps)
        for i, contact_map in enumerate(cmaps):
            contact_map[~masks[i]] = 0.

    for i in range(len(cmaps)):
        axes[i].imshow(np.array(cmaps[i]), cmap='viridis')
        axes[i].set_title(titles[i])
        axes[i].axis('off')

    if suptitle is not None:
        plt.suptitle(suptitle)
        
    plt.tight_layout()
    if save:
        plt.savefig(suptitle+'.png' if suptitle is not None else 'Comparison_default.png')
    plt.show()

def get_soft_range_matrix(contact_maps: list, mask, configs):
    stacked = np.stack(contact_maps, axis=2)  # (N, N, C)
    N = stacked.shape[0]
    C = stacked.shape[2]
    range_matrix = np.zeros((N, N, 2))

    # Extraire les limites inférieures et supérieures des configs
    lowers = np.array([cfg['lower'] for cfg in configs])
    uppers = np.array([cfg['upper'] for cfg in configs])

    for i in range(N):
        for j in range(i+1, N):
            weights = stacked[i, j]
            total = np.sum(weights)

            if total == 0:
                continue

            # Normaliser les poids (softmax facultatif)
            weights = weights / total

            # Barycentre pondéré des bornes
            lower = np.sum(weights * lowers)
            upper = np.sum(weights * uppers)

            range_matrix[i, j] = [lower, upper]
            range_matrix[j, i] = [lower, upper]

    return range_matrix

def get_range_matrix(contact_maps:list, mask, configs):

    stacked = np.stack(contact_maps, axis=2)
    range_matrix = np.zeros((contact_maps[0].shape[0], contact_maps[0].shape[1], 2))
    
    for i in range(contact_maps[0].shape[0]):
        for ii in range(contact_maps[0].shape[1]):
            if np.sum(stacked[i, ii])==0: 
                continue
            else:
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

def compute_pos_weights(contact_maps, local_mask):
    pos_weights = []
    for cmap in contact_maps:
        ratio = max(1.0, local_mask.sum() / (cmap.sum() + 1e-8))
        pos_weight = max(1.0, ratio - 1.0)
        pos_weights.append(pos_weight)
    return pos_weights

def make_local_mask(N, window):
    mask = torch.zeros((N, N), dtype=torch.bool)
    for i in range(N):
        start, end = max(0, i - window), min(N, i + window + 1)
        mask[i, start:end] = True
    return mask

def compute_logits(dists, config, sharpness):
    if config['lower'] == 0:
        return sharpness * (config['upper'] - dists) * 0.1 * torch.rand_like(dists)
    else:
        upper_logits = sharpness * (config['upper'] - dists)
        lower_logits = sharpness * (dists - config['lower'])
        return torch.min(upper_logits, lower_logits)

def apply_weighted_bce(logits, target_map, pos_weight, bce_loss_fn):
    weight_matrix = pos_weight * target_map + (1 - target_map)
    map_loss = bce_loss_fn(logits, target_map)
    return (map_loss * weight_matrix).mean()

def repulsion_term(dists, min_distance, mask):
    min_dist_violation = torch.relu(min_distance - dists)
    return (min_dist_violation[mask] ** 2).sum()

def optim_reconstruct_coords_maps(
    contact_maps, configs, hard_mask=None, sharpness=10.0, window=25,
    dim=3, lr=1e-3, max_iter=1000, prints=10, verbose=True
):
    """
    Optimizes 3D coordinates from a stack of contact maps using differentiable loss.

    Args:
        contact_maps: numpy array of shape (N, N, x)
        configs: list of configuration dicts (len must be x)
        hard_mask: optional binary mask (N, N) for positive weight selection
        sharpness: controls sigmoid steepness in loss
        window: neighborhood window for local mask if hard_mask is None
        dim: output dimension (typically 3)
        lr: learning rate
        max_iter: optimization steps
        prints: how often to log loss
        verbose: whether to print progress

    Returns:
        best_coords: optimized coordinate array (N, dim)
        best_loss: final loss
    """
    assert contact_maps.ndim == 3, "Input contact_maps must be shape (N, N, x)"
    N, _, x = contact_maps.shape
    assert len(configs) == x, "configs length must match number of maps"

    # Initialization
    coords = advanced_initialize_coords(N, dim, [contact_maps[:, :, i] for i in range(x)], configs)
    coords.requires_grad = True
    optimizer = torch.optim.Adam([coords], lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', factor=0.5, patience=50)

    contact_maps_tensor = torch.tensor(contact_maps, dtype=torch.float32)
    repulsion_mask = ~torch.eye(N, dtype=torch.bool)

    local_mask = torch.tensor(hard_mask.astype(bool)) if hard_mask is not None else make_local_mask(N, window)
    pos_weights = compute_pos_weights([contact_maps[:, :, i] for i in range(x)], local_mask)

    bce_loss_fn = torch.nn.BCEWithLogitsLoss(reduction='none')
    best_loss = float('inf')
    best_coords = coords.clone().detach()
    patience_counter = 0

    for step in range(max_iter):
        optimizer.zero_grad()
        dists = torch.cdist(coords, coords)

        # Regularization losses
        sequential_dists = torch.diag(dists, diagonal=1)
        connectivity_loss = ((sequential_dists - 3.8) ** 2).mean()
        repulsion_loss = repulsion_term(dists, min_distance=2.3, mask=repulsion_mask)

        # Main loss from all contact maps
        total_loss = torch.tensor(0.0, requires_grad=True)
        partial_losses = []
        for i in range(x):
            logits = compute_logits(dists, configs[i], sharpness)
            loss = apply_weighted_bce(logits, contact_maps_tensor[:, :, i], pos_weights[i], bce_loss_fn)
            total_loss = total_loss + loss
            partial_losses.append(loss)

        total_loss += 0.15 * repulsion_loss + 0.05 * connectivity_loss
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_([coords], 1.0)
        optimizer.step()
        scheduler.step(total_loss)

        if verbose and step % max(1, int(max_iter // prints)) == 0:
            print(f"Step {step:6d} | Total loss: {total_loss.item():8.3f} | Parts: {[f'{l.item():.2f}' for l in partial_losses]}")

        if total_loss.item() < best_loss:
            best_loss = total_loss.item()
            best_coords = coords.clone().detach()
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter > 100:
                if verbose:
                    print("Early stopping.")
                break

    return best_coords.detach().cpu().numpy(), best_loss
