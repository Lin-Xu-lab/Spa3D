import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import numba

@numba.njit("f4(f4[:], f4[:])")
def euclid_dist(t1,t2):
	sum=0
	for i in range(t1.shape[0]):
		sum+=(t1[i]-t2[i])**2
	return np.sqrt(sum)

@numba.njit("f4[:,:](f4[:,:])", parallel=True, nogil=True)
def pairwise_distance(X):
	n=X.shape[0]
	adj=np.empty((n, n), dtype=np.float32)
	for i in numba.prange(n):
		for j in numba.prange(n):
			adj[i][j]=euclid_dist(X[i], X[j])
	return adj

def calculate_adj_matrix_3D(x, y, m, x_pixel=None, y_pixel=None, m_array=None, image=None, beta=49, alpha=1, histology=False):
	#x,y,x_pixel, y_pixel are lists
	if histology:
		assert (x_pixel is not None) & (x_pixel is not None) & (image is not None)
		assert (len(x)==len(x_pixel)) & (len(y)==len(y_pixel))
		print("Calculateing adj matrix using histology image...")
		#beta to control the range of neighbourhood when calculate grey vale for one spot
		#alpha to control the color scale
		beta_half=round(beta/2)
		g=[]
		for i in range(len(x_pixel)):
			max_x=image.shape[1]
			max_y=image.shape[2]
			nbs=image[m_array[i], max(0,x_pixel[i]-beta_half):min(max_x,x_pixel[i]+beta_half+1),max(0,y_pixel[i]-beta_half):min(max_y,y_pixel[i]+beta_half+1)]
			g.append(np.mean(np.mean(nbs,axis=0),axis=0))
		c0, c1, c2=[], [], []
		for i in g:
			c0.append(i[0])
			c1.append(i[1])
			c2.append(i[2])
		c0=np.array(c0)
		c1=np.array(c1)
		c2=np.array(c2)
		print("Var of c0,c1,c2 = ", np.var(c0),np.var(c1),np.var(c2))
		c3=(c0*np.var(c0)+c1*np.var(c1)+c2*np.var(c2))/(np.var(c0)+np.var(c1)+np.var(c2))
		c4=(c3-np.mean(c3))/np.std(c3)
		z_scale=np.max([np.std(x), np.std(y)])*alpha
		z=c4*z_scale
		z=z.tolist()
		print("Var of x,y,z = ", np.var(x),np.var(y),np.var(m),np.var(z))
		X=np.array([x, y, m, z]).T.astype(np.float32)
	else:
		print("Calculateing adj matrix using xyz only...")
		X=np.array([x, y, m]).T.astype(np.float32)
		
	return pairwise_distance(X)
