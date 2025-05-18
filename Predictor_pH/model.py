import os
import pickle
import pandas as pd
import numpy as np
import torch
from torch_geometric.data import Data, Dataset
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool, GATConv

class GNN(torch.nn.Module):
    def __init__(self, num_node_features, num_edge_features, hidden_channels=64):
        super(GNN, self).__init__()
        
        # Graph convolutional layers
        self.conv1 = GATConv(num_node_features, hidden_channels)
        self.conv2 = GATConv(hidden_channels, hidden_channels)
        self.conv3 = GATConv(hidden_channels, hidden_channels)
        
        # MLP for final prediction
        self.lin1 = torch.nn.Linear(hidden_channels, hidden_channels)
        self.lin2 = torch.nn.Linear(hidden_channels, 1)
        
        # Batch normalization
        self.bn1 = torch.nn.BatchNorm1d(hidden_channels)
        self.bn2 = torch.nn.BatchNorm1d(hidden_channels)
        self.bn3 = torch.nn.BatchNorm1d(hidden_channels)
        
    def forward(self, x, edge_index, edge_attr, batch):
        # Graph convolution layers
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.bn1(x)
        x = F.dropout(x, p=0.2, training=self.training)
        
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        x = self.bn2(x)
        x = F.dropout(x, p=0.2, training=self.training)
        
        x = self.conv3(x, edge_index)
        x = F.relu(x)
        x = self.bn3(x)
        
        # Global pooling
        x = global_mean_pool(x, batch)
        
        # MLP for prediction
        x = self.lin1(x)
        x = F.relu(x)
        x = F.dropout(x, p=0.2, training=self.training)
        x = self.lin2(x)
        
        return x