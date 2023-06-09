#cd Desktop/GCN/pygcn/env2/
import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from anndata import AnnData
import torch
from sklearn.decomposition import PCA
import math
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from models import *

class Spa3D(object):
    def __init__(self):
        super(Spa3D, self).__init__()
        self.l=None

    def set_l(self, l):
        self.l=l

    def train(self,adata,adj, 
            num_pcs=50, 
            lr=0.005,
            max_epochs=2000,
            weight_decay=0,
            opt="admin",
            init_spa=True,
            init="louvain", #louvain or kmeans
            num_neighbors=10, #for louvain
            num_clusters=None, #for kmeans
            res=0.4, #for louvain
            tol=1e-3):
        self.num_pcs=num_pcs
        self.res=res
        self.lr=lr
        self.max_epochs=max_epochs
        self.weight_decay=weight_decay
        self.opt=opt
        self.init_spa=init_spa
        self.init=init
        self.num_neighbors=num_neighbors
        self.num_clusters=num_clusters
        self.res=res
        self.tol=tol
        assert adata.shape[0]==adj.shape[0]==adj.shape[1]
        pca = PCA(n_components=self.num_pcs)
        if issparse(adata.X):
            pca.fit(adata.X.A)
            embed=pca.transform(adata.X.A)
        else:
            pca.fit(adata.X)
            embed=pca.transform(adata.X)
        ###------------------------------------------###
        if self.l is None:
            raise ValueError('l should not be set before fitting the model!')
        adj_exp=np.exp(-1*(adj**2)/(2*(self.l**2)))
        #----------Train model----------
        self.model=simple_GC_DEC(embed.shape[1],embed.shape[1])
        self.model.fit(embed,adj_exp,lr=self.lr,max_epochs=self.max_epochs,weight_decay=self.weight_decay,opt=self.opt,init_spa=self.init_spa,init=self.init,num_neighbors=self.num_neighbors,num_clusters=self.num_clusters,res=self.res, tol=self.tol)
        self.embed=embed
        self.adj_exp=adj_exp

    def predict(self):
        z,q=self.model.predict(self.embed,self.adj_exp)
        y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
        # Max probability plot
        prob=q.detach().numpy()
        X=z.detach().numpy()
        return y_pred, prob, X

