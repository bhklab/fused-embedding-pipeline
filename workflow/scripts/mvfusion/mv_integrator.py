"""
Code to expand on integrao. Adapted from the original repo by James Bannon. 


original code:  Shihao Ma at WangLab, Dec. 14, 2020


"""

from typing import List, Tuple, Union, Optional, Dict
import pandas as pd
import numpy as np 
from collections import defaultdict
import snf
from .utils import _scaling_normalized_pd, _stable_normalized_pd,dist2, _stable_normalized
import time
import sys
from .models import IntegrAO
from sklearn.utils.validation import (
    check_array,
    check_symmetric,
    check_consistent_length,
)

from .training_utils import tsne_p_deep

import torch
import torch_geometric.transforms as T
import torch.nn.functional as F

class MultiViewIntegrator:
    def __init__(
            self,
            views: List[pd.DataFrame],
            view_names:List[str] = None,
            metrics: List[str] = None,
            neighborhood_size:int=None,
            mu: float = 0.5,
            norm_factor: float = 1.0,
            random_seed:int = 42,
            num_fusion_iters:int = 20,
            hidden_dim:int = 128,
            emb_dim: int = 50,
            alignment_epochs:int = 200,
            seed:int = 42
    ):
        

        """
        Initializes an object that integrates mutliple 'views' of a kind of objects
        


        Parameters
        ----------

        views : List[pd.Dataframe]
                    List of views (datasets) for consideration. 
                    
                    These are expected to be pandas dataframes where each row corresponds to a sample
                    and each column to a feature. The datasets can have different numbers of samples 
                    as well as different feature dimensions. Samples in a single view should have the
                    same number of features. 

        view_names: List[str], default None
                    List of the names for each view. These are optional but useful.


        metrics: List[str], default None
                    List of metrics to use when computing affinity for each view. Default is 'None' which uses
                    the square Euclidean distance for each. 

        neighborhood_size: int, default None
                    The cardinality of the neighborhood each node in the graph. The K in k nearest enighbors

        
        Returns
        --------
        None
        """
        
        assert len(views)>=1, "must have at least one view"
        
        self.views = views
        self.view_names = view_names
        self.mu = mu
        self.num_fusion_iters = num_fusion_iters
        self.hidden_dim = hidden_dim
        self.emb_dim = emb_dim
        self.alignment_epochs = alignment_epochs
        self.random_seed = seed
        if metrics is None:
            self.metrics = ['sqeuclidean']*len(self.views)
        else:
            assert len(metrics) == len(views), "Number of metrics and number of views must be the same or metrics must be None"
            self.metrics = metrics

        if neighborhood_size == None:
            self.neighborhood_size = int(self.views[0].shape[0] / 6)
        else:
            self.neighborhood_size = neighborhood_size
        
        self._index_views()
        self._make_affinities()
        # for k in self.affinity_matrices:
        #     print(k)
        
        self.network_fusion()
        # for k in self.fused_networks:
        #     # print(k)

        torch.manual_seed(self.random_seed)
        # sys.exit()

    
   


    def _make_affinities(self)-> None:
        self.affinity_matrices = []
       
        for i in range(0, len(self.views)):
            # print(i)
            view = self.views[i]
            metric = self.metrics[i]
            # print(self.view_names[i])
            
            # dist_mat = dist2(view.values, view.values)
            # print(self.neighborhood_size)
            if metric == 'precomputed':
                 S_mat = snf.compute.affinity_matrix(
                    view, K=self.neighborhood_size, mu=self.mu
                )
                
            else:
                # print(f"computing with metric {metric}")
                S_mat = snf.compute.make_affinity(
                    view,metric =metric, K=self.neighborhood_size, mu=self.mu
                )
               

            
            S_df = pd.DataFrame(
                data=S_mat, index=self.original_order[i], columns=self.original_order[i]
            )

            self.affinity_matrices.append(S_df)
        # for mat in self.affinity_matrices:
        #     mat.set_flags(write=1)
    

    
    def _index_views(self,verbose:bool=True)-> None:
        
        """
        Index the views for mutli-view fusion.
        """

        original_order = [0] * (len(self.views))
        dict_original_order = {}
        dict_commonSample = {}
        dict_uniqueSample = {}
        dict_commonSampleIndex = {}
        dict_sampleToIndexs = defaultdict(list)

        for i in range(0, len(self.views)):
            original_order[i] = list(self.views[i].index)
            dict_original_order[i] = original_order[i]
            for sample in original_order[i]:
                dict_sampleToIndexs[sample].append(
                    (i, np.argwhere(self.views[i].index == sample).squeeze().tolist())
                )

        for i in range(0, len(original_order)):
            for j in range(i + 1, len(original_order)):
                commonList = list(set(original_order[i]).intersection(original_order[j]))
                if verbose:
                    print("Common sample between view{} and view{}: {}".format(i, j, len(commonList)))
                dict_commonSample.update(
                    dict_commonSample.fromkeys([(i, j), (j, i)], commonList)
                )
                dict_commonSampleIndex[(i, j)] = [
                    np.argwhere(self.views[i].index == x).squeeze().tolist()
                    for x in commonList
                ]
                dict_commonSampleIndex[(j, i)] = [
                    np.argwhere(self.views[j].index == x).squeeze().tolist()
                    for x in commonList
                ]

                dict_uniqueSample[(i, j)] = list(
                    set(original_order[i]).symmetric_difference(commonList)
                )
                dict_uniqueSample[(j, i)] = list(
                    set(original_order[j]).symmetric_difference(commonList)
                )

        self.dicts_common = dict_commonSample
        self.dicts_commonIndex = dict_commonSampleIndex
        self.dict_sampleToIndexs = dict_sampleToIndexs
        self.dicts_unique = dict_uniqueSample
        self.original_order = original_order
        self.dict_original_order = dict_original_order


    def network_fusion(self,verbose:bool = True)->None:
        """
        Performs Patient Graph Fusion on pairwise affinity matrices

        """
        if verbose:
            print("Start applying diffusion!")

        start_time = time.time()
        self.fused_networks = []
        self.fused_affinities = [0] * len(self.affinity_matrices)
        
        
        # First, normalize different networks to avoid scale problems, it is compatible with pandas dataframe
        for n, mat in enumerate(self.affinity_matrices):
            # print(n)
            # print(mat)
            # sys.exit()
        
            
            # normalize affinity matrix based on strength of edges
            # mat = mat / np.nansum(mat, axis=1, keepdims=True) 
            # chunk =  _stable_normalized_pd(mat.copy())

            self.affinity_matrices[n] = _stable_normalized_pd(mat)
            # print("Self affinities")
           
            # aff[n] = check_symmetric(mat, raise_warning=False)

            # apply KNN threshold to normalized affinity matrix
            # We need to crop the intersecting samples from newW matrices
            neighborhood_size = min(int(self.neighborhood_size), mat.shape[0])
            self.fused_affinities[n] = snf.compute._find_dominate_set(self.affinity_matrices[n], neighborhood_size)
           

        # If there is only one view, return it
        # print(self.fused_affinities)
        
        if len(self.affinity_matrices) == 1:
            if verbose:
                print("Only one view, return it directly")
            self.fused_networks = self.fused_affinities
            # print(self.fused_networks)
            

        else:
            for iteration in range(self.num_fusion_iters):

                # Make a copy of the aff matrix for this iteration
                # goal is to update aff[n], but it is the average of all the defused matrices
                # Make a copy of aff[n], and set it to 0
                aff_next = []
                for k in range(len(self.affinity_matrices)):
                    aff_temp = self.affinity_matrices[k].copy()
                    aff_temp = pd.DataFrame(
                        np.zeros(shape=self.affinity_matrices[k].shape),
                        columns = self.affinity_matrices[k].columns,
                        index=self.affinity_matrices[k].index)
                    
                    # for col in aff_temp.columns:
                    #     aff_temp[col].values[:] = 0
                    aff_next.append(aff_temp)

                for n, mat in enumerate(self.affinity_matrices):
                    nzW = self.fused_affinities[n]  
                
                    

                    for j, mat_tofuse in enumerate(self.affinity_matrices):
                        if n == j:
                            continue

                        raw_fusion_matrix =  pd.DataFrame(
                            data=np.identity(nzW.shape[0]),
                            index=self.original_order[n],
                            columns=self.original_order[n],
                        ) + mat_tofuse.loc[sorted(self.dicts_common[(j, n)]),sorted(self.dicts_common[(j, n)])].copy()
                        
                        raw_fusion_matrix.fillna(0.0, inplace=True)
                        
                        fusion_matrix = _scaling_normalized_pd(raw_fusion_matrix, ratio=len(self.dicts_common[(j, n)])/raw_fusion_matrix.shape[0])
                    
                        fusion_matrix = check_symmetric(fusion_matrix, raise_warning=False)
                        fusion_matrix = fusion_matrix.reindex(self.original_order[n], axis=1)
                        fusion_matrix = fusion_matrix.reindex(self.original_order[n], axis=0)

                        # Now we are ready to do the diffusion
                        nzW_T = np.transpose(nzW)
                        aff0_temp = nzW.dot(
                            fusion_matrix.dot(nzW_T)
                        )  # Matmul is not working, but .dot() is good

                        aff0_temp = _stable_normalized_pd(aff0_temp)
                        
                        

                        aff_next[n] = np.add(aff0_temp, aff_next[n])

                    aff_next[n] = np.divide(aff_next[n], len(self.affinity_matrices) - 1)
                    # aff_next[n] = _stable_normalized_pd(aff_next[n])

                # put the value in aff_next back to aff
                for k in range(len(self.affinity_matrices)):
                    self.affinity_matrices[k] = aff_next[k]

            for n, mat in enumerate(self.affinity_matrices):
                self.affinity_matrices[n] = _stable_normalized_pd(mat.copy())
            
            self.fused_networks = self.affinity_matrices
            end_time = time.time()
            if verbose:
                print("Diffusion ends! Times: {}s".format(end_time - start_time))
        


    def neural_integration(self)->None:
    
        datasets_val = [x.values for x in self.views]
        fused_networks_val = [x.values.copy() for x in self.fused_networks]
        # print(fused_networks_val)
        S_final, self.models = tsne_p_deep(
            self.dicts_commonIndex,
            self.dict_sampleToIndexs,
            datasets_val,
            P=fused_networks_val,
            neighbor_size=self.neighborhood_size,
            embedding_dims=self.emb_dim,
            alighment_epochs=self.alignment_epochs,
            seed = self.random_seed
        )

        self.final_embeds = pd.DataFrame(
            data=S_final, index=self.dict_sampleToIndexs.keys()
        )
        self.final_embeds.sort_index(inplace=True)

        # calculate the final similarity graph
        dist_final = dist2(self.final_embeds.values, self.final_embeds.values)
        Wall_final = snf.compute.affinity_matrix(
            dist_final, K=self.neighborhood_size, mu=self.mu
        )

        Wall_final = _stable_normalized(Wall_final)

        return self.final_embeds, Wall_final, self.models
        

    def get_fused_matrices(self):
        fused_mats = {}
        for i in range(len(self.fused_affinities)):
            mat = self.fused_affinities[i]
            idx = self.original_order[i]
            view = self.view_names[i]
            labeled_mat = pd.DataFrame(mat,index=idx,columns=idx)
            fused_mats[view] = labeled_mat
        return fused_mats


    def get_affinity_matrices(self):
        aff_mats = {}
        for i in range(len(self.affinity_matrices)):
            mat = self.affinity_matrices[i]
            idx = self.original_order[i]
            view = self.view_names[i]
            labeled_mat = pd.DataFrame(mat,index=idx,columns=idx)
            aff_mats[view] = labeled_mat
        return aff_mats