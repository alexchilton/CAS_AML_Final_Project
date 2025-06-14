import plotly.graph_objects as go
from plotly.subplots import make_subplots
from matplotlib import cm
import numpy as np
import torch
import itertools

def title(text):
    print (f'\033[1;31m{text}\033[0m')

def rotate_and_mirror(structure):

    theta = np.pi / 1.5
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [0, 0, 1]
    ])
    
    mirror_axes = (True, False, True)
    reflection = np.diag([ -1 if mirror else 1 for mirror in mirror_axes ])
    mirrored = structure @ reflection
    rotated = mirrored @ rotation_matrix.T  # Note: transpose because coords are row vectors
    return rotated

def huber_loss(residuals, delta=1.0):
    return np.where(residuals < delta,
                    0.5 * residuals ** 2,
                    delta * (residuals - 0.5 * delta))

def align_structures(A, B, ref_slice=slice(0, 3)):

    assert A.shape[1] == 3 and B.shape[1] == 3, "Structures must be Nx3"
    assert A.shape[0] == B.shape[0], "Structures must have same number of points"
    
    mirror_options = list(itertools.product([1, -1], repeat=3))
    best_rmsd = float('inf')
    best_result = None
    
    A_ref = A[ref_slice]
    for mirror in mirror_options:
        B_mirrored = B * np.array(mirror)
        B_ref = B_mirrored[ref_slice]

        # Center reference points
        A_centered = A_ref - A_ref.mean(axis=0)
        B_centered = B_ref - B_ref.mean(axis=0)

        # Kabsch alignment
        H = B_centered.T @ A_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T

        # Apply transformation to all of B
        B_rotated = (B_mirrored - B_ref.mean(axis=0)) @ R
        B_aligned = B_rotated + A_ref.mean(axis=0)

        # Compute RMSD
        rmsd = np.sqrt(np.mean(np.sum((A - B_aligned) ** 2, axis=1)))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_result = (B_aligned, R, A_ref.mean(axis=0) - B_ref.mean(axis=0) @ R, mirror, rmsd)

    return best_result  # B_aligned, R, t, mirror, rmsd
    
def plotly_2_graphs(data: list, title='Default title', labels=[]):
    
    num_points = len(data[0])
    colors = cm.rainbow(np.linspace(0, 1, num_points))
    
    # Create subplot layout (1 row, 2 columns, with 3D specs)
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]]
        # subplot_titles=("Original structure", "Reconstructed structure")
    )

    # First plot (original)
    d = data[0]
    x, y, z = d[:,0],d[:,1],d[:,2]
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers+lines',
            marker=dict(size=3, color=colors),
            line=dict(width=1, color=colors),
            name='Original'
        ),
        row=1, col=1
    )
    
    # Second plot (reconstructed or same again)
    d = data[1]
    x, y, z = d[:,0],d[:,1],d[:,2]
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers+lines',
            marker=dict(size=3, color=colors),
            line=dict(width=1, color=colors),
            name='Reconstructed'
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        scene1=dict(  # Notez "scene1" pour le premier graphique
            domain=dict(x=[0.0, 0.48], y=[0.0,0.9]),
            xaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='grey',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            yaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='grey',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            zaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='grey',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            bgcolor='rgba(0,0,0,0)',
            camera=dict(
                eye=dict(x=1.0, y=1.0, z=1.0),  # Augmente le zoom en réduisant les valeurs
                center=dict(x=0, y=0, z=0)
            )
        ),
        
        scene2=dict(  # Notez "scene2" pour le deuxième graphique
            domain=dict(x=[0.52, 1.], y=[0.0,0.9]),
            xaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='#999999',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            yaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='#999999',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            zaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                visible=False,           # Rendre l'axe visible
                showline=True,          # Afficher la ligne de l'axe
                linecolor='#999999',  # Couleur grise claire pour la ligne
                linewidth=1,            # Épaisseur de ligne minimale
                showspikes=False        # Pas de "spikes" lors du survol
            ),
            bgcolor='rgba(0,0,0,0)',
            camera=dict(
                eye=dict(x=1.0, y=1.0, z=1.0),  # Augmente le zoom en réduisant les valeurs
                center=dict(x=0, y=0, z=0)
            ),
        ),
        
        height=400,
        width=800,
        #title_text=title,
        showlegend=False,
        
        # Ajout d'un outline pour la figure entière
        paper_bgcolor='white',
        plot_bgcolor='white',
        
        margin=dict(
            l=5,   # marge gauche
            r=5,   # marge droite
            t=50,  # marge supérieure (pour le titre)
            b=5    # marge inférieure
        ),

        font=dict(
            family="Menlo, monospace",
            size=10,        
            color="#666666"   
        ),

        # Ajouter une bordure autour de la figure
        shapes=[
            dict(
                type="rect",
                xref="paper",
                yref="paper",
                x0=0,
                y0=0,
                x1=0.48,
                y1=0.9,
                line=dict(
                    color="#999999",
                    width=1,
                )
            ),
            
            dict(
                type="rect",
                xref="paper",
                yref="paper",
                x0=0.52,
                y0=0,
                x1=1,
                y1=0.9,
                line=dict(
                    color="#999999",
                    width=1,
                )
            )  
        ],

        annotations=[
            dict(
                text=title,
                x=0,
                y=1.15,
                xref='paper',
                yref='paper',
                showarrow=False,
                xanchor='left',
                align='left',
                font=dict(size=16, family="Menlo, monospace")
            ),
            dict(
                text="Original structure",
                x=0.25,  # Center of subplot 1
                y=0.05,   # Below the plot
                xref="paper",
                yref="paper",
                xanchor="center",
                showarrow=False,
                font=dict(size=14, family="Menlo, monospace")
            ),
            dict(
                text="Recovered structure",
                x=0.75,  # Center of subplot 2
                y=0.05,
                xref="paper",
                yref="paper",
                xanchor="center",
                showarrow=False,
                font=dict(size=14, family="Menlo, monospace")
            )
        ] + [
            dict(
                text='⚿ '+text,
                x=0.,
                y=1.08 - i * 0.05,
                xref='paper',
                yref='paper',
                showarrow=False,
                xanchor='left',
                align='left',
                font=dict(size=12, family="Menlo, monospace")
            ) for i, text in enumerate(labels)
        ]
    )
    
    
    fig.show()