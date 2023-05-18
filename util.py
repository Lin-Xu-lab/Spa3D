import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import os
from anndata import AnnData,read_csv,read_text,read_mtx
from scipy.sparse import issparse
import random
import torch
from Spa3D import Spa3D

def preprocessing_filter_by_gene_number(adata,min_counts=None,max_counts=None,min_cells=10,max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_cells=min_cells)[0]) if min_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_cells=max_cells)[0]) if max_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_var(id_tmp)

def prefilter_specialgenes(adata,Gene1Pattern="ERCC",Gene2Pattern="MT-"):
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)

def calculate_p(adj, l):
    adj_exp=np.exp(-1*(adj**2)/(2*(l**2)))
    return np.mean(np.sum(adj_exp,1))-1

def search_decay_value_3D(adj, p = 0.5, start = 0.01, end = 1000, tol = 0.000001, max_run = 100):
    run=0
    p_low=calculate_p(adj, start)
    p_high=calculate_p(adj, end)
    if p_low>p+tol:
        print("l not found, try smaller start point.")
        return None
    elif p_high<p-tol:
        print("l not found, try bigger end point.")
        return None
    elif  np.abs(p_low-p) <=tol:
        print("We recommend to use decay_value = ", str(start))
        return start
    elif  np.abs(p_high-p) <=tol:
        print("We recommend to use decay_value = ", str(end))
        return end
    while (p_low+tol)<p<(p_high-tol):
        run+=1
        print("Iteration ID: "+str(run)+": decay_value ["+str(start)+", "+str(end)+"], p ["+str(p_low)+", "+str(p_high)+"]")
        if run >max_run:
            print("Exact l not found, closest values are:\n"+"decay_value = "+str(start)+": "+"p = "+str(p_low)+"\nl="+str(end)+": "+"p = "+str(p_high))
            return None
        mid=(start+end)/2
        p_mid=calculate_p(adj, mid)
        if np.abs(p_mid-p)<=tol:
            print("We recommend to use decay_value = ", str(mid))
            return mid
        if p_mid<=p:
            start=mid
            p_low=p_mid
        else:
            end=mid
            p_high=p_mid

def search_proper_resolution_3D(adata, adj, decay_value, num_clusters, num_neighbors = 20, start = 0.7, step = 0.1, tol = 5e-6, lr = 0.05, max_epochs = 20,  r_seed = 100, t_seed = 100, n_seed = 100, max_run = 10):
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    res=start
    print("Start at res = ", res, "step = ", step)
    clf=Spa3D()
    clf.set_l(decay_value)
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=tol, lr=lr, max_epochs=max_epochs, num_neighbors = num_neighbors)
    y_pred, _, _=clf.predict()
    old_num=len(set(y_pred))
    print("Res = ", res, "Num of clusters = ", old_num)
    run=0
    while old_num!=num_clusters:
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        old_sign=1 if (old_num<num_clusters) else -1
        clf=Spa3D()
        clf.set_l(decay_value)
        clf.train(adata,adj,init_spa=True,init="louvain",res=res+step*old_sign, tol=tol, lr=lr, max_epochs=max_epochs, num_neighbors = num_neighbors)
        y_pred, _, _=clf.predict()
        new_num=len(set(y_pred))
        print("Res = ", res+step*old_sign, "Num of clusters = ", new_num)
        if new_num==num_clusters:
            res=res+step*old_sign
            print("recommended res = ", str(res))
            return res
        new_sign=1 if (new_num<num_clusters) else -1
        if new_sign==old_sign:
            res=res+step*old_sign
            print("Res changed to", res)
            old_num=new_num
        else:
            step=step/2
            print("Step changed to", step)
        if run >max_run:
            print("Exact resolution not found")
            print("Recommended resolution = ", str(res))
            return res
        run+=1
    print("recommended resolution = ", str(res))
    return res


def refine(sample_id, pred, dis, shape="hexagon"):
    refined_pred=[]
    pred=pd.DataFrame({"pred": pred}, index=sample_id)
    dis_df=pd.DataFrame(dis, index=sample_id, columns=sample_id)
    if shape=="hexagon":
        num_nbs=6 
    elif shape=="square":
        num_nbs=4
    elif shape=="eight":
    	num_nbs=8
    elif shape=="ten":
    	num_nbs=10
    elif shape=="12":
    	num_nbs=12
    elif shape=="14":
    	num_nbs=14
    elif shape=="16":
    	num_nbs=16
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=dis_df.loc[index, :].sort_values()
        nbs=dis_tmp[0:num_nbs+1]
        nbs_pred=pred.loc[nbs.index, "pred"]
        self_pred=pred.loc[index, "pred"]
        v_c=nbs_pred.value_counts()
        if (v_c.loc[self_pred]<num_nbs/2) and (np.max(v_c)>num_nbs/2):
            refined_pred.append(v_c.idxmax())
        else:           
            refined_pred.append(self_pred)
    return refined_pred
