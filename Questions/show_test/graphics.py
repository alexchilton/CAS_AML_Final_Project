import plotly.graph_objects as go
from plotly.subplots import make_subplots
from matplotlib import cm
import numpy as np

def plotly_2_graphs(data: list, title='Default title'):
    
    num_points = len(data[0])
    colors = cm.rainbow(np.linspace(0, 1, num_points))
    
    # Create subplot layout (1 row, 2 columns, with 3D specs)
    fig = make_subplots(
        rows=1, cols=2,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}]],
        subplot_titles=("Original", "Reconstructed")
    )

    # First plot (original)
    d = data[0]
    x, y, z = d[:,0],d[:,1],d[:,2]
    fig.add_trace(
        go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers+lines',
            marker=dict(size=4, color=colors),
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
            marker=dict(size=4, color=colors),
            line=dict(width=1, color=colors),
            name='Reconstructed'
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        scene1=dict(  # Notez "scene1" pour le premier graphique
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
        height=500,
        width=1000,
        title_text=title,
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
        # Ajouter une bordure autour de la figure
        shapes=[
            dict(
                type="rect",
                xref="paper",
                yref="paper",
                x0=0,
                y0=0,
                x1=1,
                y1=1,
                line=dict(
                    color="grey",
                    width=1,
                )
            )
        ]
    )
    
    
    fig.show()