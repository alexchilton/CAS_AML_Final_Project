import os
import pickle
import pandas as pd
import numpy as np
import torch
from torch_geometric.data import Data, Dataset
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool, GATConv
from model import GNN


class EnsembleGNN(torch.nn.Module):
    def __init__(self, num_node_features, num_edge_features, num_chemical_features, hidden_channels=64):
        super(EnsembleGNN, self).__init__()
        
        # Graph neural network component
        self.gnn = GNN(num_node_features, num_edge_features, hidden_channels)
        
        # Chemical features processing
        self.chemical_nn = torch.nn.Sequential(
            torch.nn.Linear(num_chemical_features, hidden_channels),
            torch.nn.ReLU(),
            torch.nn.BatchNorm1d(hidden_channels),
            torch.nn.Dropout(0.2),
            torch.nn.Linear(hidden_channels, hidden_channels)
        )
        
        # Fusion and final prediction
        self.fusion = torch.nn.Sequential(
            torch.nn.Linear(hidden_channels * 2, hidden_channels),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.2),
            torch.nn.Linear(hidden_channels, 1)
        )
        
    def forward(self, x, edge_index, edge_attr, batch, chemical_features):
        # Process graph structure
        gnn_out = self.gnn(x, edge_index, edge_attr, batch)
        
        # Process chemical features
        chem_out = self.chemical_nn(chemical_features)
        
        # Combine features
        combined = torch.cat([gnn_out, chem_out], dim=1)
        
        # Make final prediction
        out = self.fusion(combined)
        
        return out