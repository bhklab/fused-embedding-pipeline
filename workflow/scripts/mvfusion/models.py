from typing import List, Tuple, Union, Optional, Dict

import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import numpy as np
from torch_geometric.nn import GraphSAGE




class IntegrAO(nn.Module):
    def __init__(self, 
                in_channels:List[int],
                hidden_channels:int = 128, 
                out_channels:int = 50,
                num_layers:int =2):
        
        """
        Initializes an object that integrates mutliple 'views' of a kind of objects
        


        Parameters
        ----------

        in_channels : List[int]
                   A list of feature dimensions for each view
                    
                    These are expected to be pandas dataframes where each row corresponds to a sample
                    and each column to a feature. The datasets can have different numbers of samples 
                    as well as different feature dimensions. Samples in a single view should have the
                    same number of features. 

        hidden_channels: int, default 128
                    The hidden dimension for each GNN layer


        out_channels: int, default 50
                    Final embedding dimension

        num_layers: int, default 2
                    The number of hidden layers per gNN

        """
        super(IntegrAO, self).__init__()
        
        self.in_channels = in_channels 
        self.hidden_channels = hidden_channels
        self.output_dim = out_channels
        self.num_layers = num_layers

        num = len(in_channels)
        
        
        view_networks = []

        for i in range(len(self.in_channels)):
            gnn= GraphSAGE(
                in_channels=self.in_channels[i], 
                hidden_channels=self.hidden_channels, 
                num_layers=self.num_layers, 
                out_channels=self.output_dim,
                project=False)
            
            view_networks.append(gnn)

        self.view_networks = nn.ModuleList(view_networks)

        self.projection_layer = nn.Sequential(
            nn.Linear(self.output_dim, self.output_dim),
            nn.BatchNorm1d(self.output_dim),
            nn.LeakyReLU(0.1, True),
            nn.Linear(self.output_dim, self.output_dim),
        )

    def forward(self, x_dict, edge_index_dict):
        z_all = {}
        for domain in x_dict.keys():
            z = self.view_networks[domain](x_dict[domain], edge_index_dict[domain])
            z = self.projection_layer(z)
            z_all[domain] = z
        return z_all
