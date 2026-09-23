from torch_geometric.data import InMemoryDataset, Data
from snf.compute import _find_dominate_set
from sklearn.utils.validation import (
    check_array,
    check_symmetric,
    check_consistent_length,
    
)
import networkx as nx
import numpy as np
import torch

# custom dataset
class GraphDataset(InMemoryDataset):
    """
    Data class for storing the views as graphs for use in PytorchGeometric
    """

    def __init__(self,
                neighbor_size, 
                feature_matrix, 
                edge_weights,
                transform=None):
        
        
        super(GraphDataset, self).__init__('.', transform, None, None)
        

       

        # compute the maximum neighborhood size
        neighbor_size = min(int(neighbor_size), edge_weights.shape[0])

        # preprocess the input into a pyg graph
        network = _find_dominate_set(edge_weights, K=neighbor_size)
        network = check_symmetric(edge_weights, raise_warning=False).copy()
        network[network > 0.0] = 1.0
        G = nx.from_numpy_array(edge_weights)

        # create edge index from 
        adj = nx.to_scipy_sparse_array(G).tocoo() 
        row = torch.from_numpy(adj.row.astype(np.int64)).to(torch.long)
        col = torch.from_numpy(adj.col.astype(np.int64)).to(torch.long)
        edge_index = torch.stack([row, col], dim=0)

        data = Data(edge_index=edge_index)
        data.num_nodes = G.number_of_nodes()
        
        # embedding 
        feature_matrix = feature_matrix.astype(np.float32)
        data.x = torch.from_numpy(feature_matrix.copy()).type(torch.float32)
        
        self.data, self.slices = self.collate([data])

    def _download(self):
        return

    def _process(self):
        return

    def __repr__(self):
        return '{}()'.format(self.__class__.__name__)
   
