import sys

import torch
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Dict, Any, Optional, Tuple
import os
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from torch_geometric.data import Data
import networkx as nx
# Add the preprocessing directory to the path
sys.path.append(os.path.abspath('../model'))
from GATVAE import GATVAE


class GenerationVisualizer:
    """
    Class for visualizing generated protein structures and latent space representations.
    """

    def __init__(self, model: GATVAE, device='cpu', save_dir='./generation_visualizations'):
        """
        Initialize the generation visualizer.

        Args:
            model: Trained GATVAE model
            device: Device to run computations on
            save_dir: Directory to save visualizations
        """
        self.model = model
        self.device = device
        self.save_dir = save_dir

        # Create directory if it doesn't exist
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        # Move model to device and set to evaluation mode
        self.model.to(device)
        self.model.eval()

    def plot_protein_graph(self, data: Data, title: str = "Protein Graph",
                           show: bool = True, save: bool = False,
                           filename: Optional[str] = None) -> nx.Graph:
        """
        Visualize a protein graph structure in 3D.

        Args:
            data: Protein graph data
            title: Plot title
            show: Whether to display the plot
            save: Whether to save the plot to disk
            filename: Filename for saved plot (if None, use title)

        Returns:
            NetworkX graph object
        """
        # Convert PyG data to NetworkX for visualization
        G = self._pyg_to_networkx(data)

        # Extract 3D coordinates from the graph
        pos = {}
        # Extract coordinates from node features (assuming they're at standard position)
        coord_start_idx = 21 + 7  # After amino acid and secondary structure features

        # Debug: print shape of node features tensor
        print(f"Node features tensor shape: {data.x.shape}")

        for i, node in enumerate(G.nodes()):
            # Debug: print node features at coordinate positions
            if i < 5:  # Print just a few nodes to avoid flooding output
                print(f"Node {i} coordinates: {data.x[i, coord_start_idx:coord_start_idx+3].cpu().numpy()}")
            # Get coordinates - assuming they're stored in the node features
            coords = data.x[i, coord_start_idx:coord_start_idx+3].cpu().numpy()
            pos[node] = coords

        # Debug: print the first few positions  # Debug: check if all coordinates are unique
        unique_coords = set(tuple(p) for p in pos.values())
        print(f"Number of unique coordinate positions: {len(unique_coords)} out of {len(pos)} nodes")
        # Create 3D plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Draw nodes
        for node, (x, y, z) in pos.items():
            ax.scatter(x, y, z, c='blue', s=50, alpha=0.7)

        # Draw edges
        for u, v in G.edges():
            x = [pos[u][0], pos[v][0]]
            y = [pos[u][1], pos[v][1]]
            z = [pos[u][2], pos[v][2]]
            ax.plot(x, y, z, c='gray', alpha=0.5)

        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Make the panes transparent
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        if save:
            if filename is None:
                filename = title.replace(" ", "_").lower() + ".png"
            plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

        return G

    def plot_interpolation_sequence(self, interpolated_graphs: List[Data],
                                    show: bool = True, save: bool = False,
                                    filename: str = "interpolation_sequence.png"):
        """
        Visualize a sequence of interpolated protein structures.

        Args:
            interpolated_graphs: List of interpolated protein graphs
            show: Whether to display the plot
            save: Whether to save the plot to disk
            filename: Filename for saved plot
        """
        n_graphs = len(interpolated_graphs)
        fig = plt.figure(figsize=(15, 5 * ((n_graphs + 2) // 3)))

        # Create subplot for each protein in the interpolation
        for i, graph in enumerate(interpolated_graphs):
            ax = fig.add_subplot(((n_graphs + 2) // 3), 3, i+1, projection='3d')

            # Use the helper method to plot the protein
            self._plot_protein_data(graph, ax, title=f"Step {i+1}")

        plt.tight_layout()

        if save:
            plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

    def _plot_protein_data(self, data: Data, ax, title: str = ""):
        """
        Helper method to plot a protein graph on a given axis.

        Args:
            data: Protein graph data
            ax: Matplotlib axis to plot on
            title: Plot title
        """
        # Extract 3D coordinates (assuming standard feature position)
        coord_start_idx = 21 + 7  # After amino acid and secondary structure features

        # Get node coordinates
        coords = data.x[:, coord_start_idx:coord_start_idx+3].cpu().numpy()

        # Plot nodes
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c='blue', s=20, alpha=0.7)

        # Plot edges
        edge_index = data.edge_index.cpu().numpy()
        for j in range(edge_index.shape[1]):
            src, dst = edge_index[0, j], edge_index[1, j]
            ax.plot([coords[src, 0], coords[dst, 0]],
                    [coords[src, 1], coords[dst, 1]],
                    [coords[src, 2], coords[dst, 2]], c='gray', alpha=0.3)

        ax.set_title(title)

        # Make the panes transparent
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        # Remove axis ticks for cleaner visualization
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    def _pyg_to_networkx(self, data: Data) -> nx.Graph:
        """
        Convert PyTorch Geometric data to NetworkX graph.

        Args:
            data: PyTorch Geometric data

        Returns:
            NetworkX graph
        """
        G = nx.Graph()

        # Add nodes
        for i in range(data.x.size(0)):
            G.add_node(i)

        # Add edges
        edge_index = data.edge_index.cpu().numpy()
        for i in range(edge_index.shape[1]):
            src, dst = edge_index[0, i], edge_index[1, i]
            G.add_edge(src, dst)

        return G  # Return the full graph, not just G.nodes
def compare_original_vs_generated(self, original_data: Data, generated_data: Data,
                                  show: bool = True, save: bool = False,
                                  filename: str = "original_vs_generated.png"):
    """
    Compare original and generated protein structures side by side.

    Args:
        original_data: Original protein graph
        generated_data: Generated protein graph
        show: Whether to display the plot
        save: Whether to save the plot to disk
        filename: Filename for saved plot
    """
    fig = plt.figure(figsize=(15, 7))

    # Plot original protein
    ax1 = fig.add_subplot(121, projection='3d')
    self._plot_protein_data(original_data, ax1, title="Original Protein")

    # Plot generated protein
    ax2 = fig.add_subplot(122, projection='3d')
    self._plot_protein_data(generated_data, ax2, title="Generated Protein")

    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def visualize_latent_space_neighborhood(self, center_z: torch.Tensor,
                                        generated_graphs: List[Data],
                                        show: bool = True, save: bool = False,
                                        filename: str = "latent_space_neighborhood.png"):
    """
    Visualize proteins generated from a neighborhood in latent space.

    Args:
        center_z: Center latent vector
        generated_graphs: List of protein graphs generated around center_z
        show: Whether to display the plot
        save: Whether to save the plot to disk
        filename: Filename for saved plot
    """
    n_graphs = len(generated_graphs)
    fig = plt.figure(figsize=(15, 5 * ((n_graphs + 3) // 4)))

    # First, plot the center protein
    ax_center = fig.add_subplot(((n_graphs + 4) // 4), 4, 1, projection='3d')

    # Assuming we have the center protein as the first element
    self._plot_protein_data(generated_graphs[0], ax_center, title="Center Protein")

    # Plot variations around the center
    for i in range(1, n_graphs):
        ax = fig.add_subplot(((n_graphs + 4) // 4), 4, i+1, projection='3d')
        self._plot_protein_data(generated_graphs[i], ax, title=f"Variation {i}")

    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def plot_latent_trajectories(self, latent_vectors: List[torch.Tensor],
                             labels: Optional[List[str]] = None,
                             show: bool = True, save: bool = False,
                             filename: str = "latent_trajectories.png"):
    """
    Visualize trajectories in the latent space.

    Args:
        latent_vectors: List of latent vectors representing the trajectory
        labels: Optional list of labels for different trajectories
        show: Whether to display the plot
        save: Whether to save the plot to disk
        filename: Filename for saved plot
    """
    # Reduce dimensionality to 2D for visualization
    all_vectors = torch.cat([v.cpu() for v in latent_vectors], dim=0)
    pca = PCA(n_components=2)
    reduced_vecs = pca.fit_transform(all_vectors.numpy())

    # Reorganize the reduced vectors back into trajectories
    trajectories = []
    offset = 0
    for vec in latent_vectors:
        trajectory_len = vec.shape[0]
        trajectories.append(reduced_vecs[offset:offset+trajectory_len])
        offset += trajectory_len

    # Plot the trajectories
    plt.figure(figsize=(10, 8))

    colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories)))

    for i, traj in enumerate(trajectories):
        label = labels[i] if labels and i < len(labels) else f"Trajectory {i+1}"
        plt.plot(traj[:, 0], traj[:, 1], 'o-', c=colors[i], label=label)

        # Mark start and end with special markers
        plt.plot(traj[0, 0], traj[0, 1], 'o', c=colors[i], markersize=10)
        plt.plot(traj[-1, 0], traj[-1, 1], 's', c=colors[i], markersize=10)

    plt.title("Latent Space Trajectories")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if save:
        plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def visualize_optimization_progress(self, latent_vectors: List[torch.Tensor],
                                    objective_values: List[float],
                                    template_data: Data,
                                    num_samples: int = 5,
                                    show: bool = True, save: bool = False,
                                    filename: str = "optimization_progress.png"):
    """
    Visualize the optimization progress with objective function values.

    Args:
        latent_vectors: List of latent vectors during optimization
        objective_values: List of objective function values
        template_data: Template graph structure to use for reconstruction
        num_samples: Number of protein samples to show
        show: Whether to display the plot
        save: Whether to save the plot to disk
        filename: Filename for saved plot
    """
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(15, 10))

    # Plot objective function values
    ax1 = fig.add_subplot(211)
    ax1.plot(objective_values, 'b-', label='Objective Value')
    ax1.set_title("Optimization Progress")
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Objective Value")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Sample iterations to visualize
    if len(latent_vectors) > num_samples:
        sample_indices = np.linspace(0, len(latent_vectors)-1, num_samples, dtype=int)
    else:
        sample_indices = range(len(latent_vectors))

    # Plot sampled proteins
    for subplot_idx, iteration_idx in enumerate(sample_indices):
        # Create 3D subplot
        ax = fig.add_subplot(2, num_samples, num_samples + subplot_idx + 1, projection='3d')

        # Generate protein from latent vector at this iteration
        with torch.no_grad():
            z = latent_vectors[iteration_idx].to(self.device)
            template_data_local = template_data.to(self.device)
            x_reconstructed = self.model.decode(z.unsqueeze(0), template_data_local)

            # Create a protein data object for visualization
            protein_data = Data(
                x=x_reconstructed,
                edge_index=template_data.edge_index,
                edge_attr=template_data.edge_attr
            )

        # Plot the protein
        self._plot_protein_data(protein_data, ax, title=f"Iteration {iteration_idx}")

    plt.tight_layout()

    if save:
        plt.savefig(os.path.join(self.save_dir, filename), dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def _plot_protein_data(self, data: Data, ax, title: str = ""):
    """
    Helper method to plot a protein graph on a given axis.

    Args:
        data: Protein graph data
        ax: Matplotlib axis to plot on
        title: Plot title
    """
    # Extract 3D coordinates (assuming standard feature position)
    coord_start_idx = 21 + 7  # After amino acid and secondary structure features

    # Get node coordinates
    coords = data.x[:, coord_start_idx:coord_start_idx+3].cpu().numpy()

    # Plot nodes
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c='blue', s=20, alpha=0.7)

    # Plot edges
    edge_index = data.edge_index.cpu().numpy()
    for j in range(edge_index.shape[1]):
        src, dst = edge_index[0, j], edge_index[1, j]
        ax.plot([coords[src, 0], coords[dst, 0]],
                [coords[src, 1], coords[dst, 1]],
                [coords[src, 2], coords[dst, 2]], c='gray', alpha=0.3)

    ax.set_title(title)

    # Make the panes transparent
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Remove axis ticks for cleaner visualization
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

def _pyg_to_networkx(self, data: Data) -> nx.Graph:
    """
    Convert PyTorch Geometric data to NetworkX graph.

    Args:
        data: PyTorch Geometric data

    Returns:
        NetworkX graph
    """
    G = nx.Graph()

    # Add nodes
    for i in range(data.x.size(0)):
        G.add_node(i)

    # Add edges
    edge_index = data.edge_index.cpu().numpy()
    for i in range(edge_index.shape[1]):
        src, dst = edge_index[0, i], edge_index[1, i]
        G.add_edge(src, dst)

    return G