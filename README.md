# 1. Introduction

Spa3D is an algorithm that utilizes spatial pattern enhancement and graph convolutional neural network model to reconstruct 3D-based spatial structure from multiple spatially resolved transcriptomics (SRT) slides with 3D spatial coordinates. 

Paper: 3D reconstruction of spatial transcriptomics with spatial pattern enhanced graph convolutional neural network

Keywords: Spatial transcriptomics, 3D reconstruction algorithm, graph convolutional network, spatial patten enhancement.

# 2. Result

Below is an example of the Spa3D for the Human embryonic heart dataset data.

![Fig](/images/Figure1-1.png)

Inserting these data into a 3D human heart model can result in:

![Fig](/images/Figure1-2.png)
    
# 3. Environment setup and code compilation

__3.1. Download the package__

The package can be downloaded by running the following command in the terminal:
```
git clone https://github.com/Lin-Xu-lab/Spa3D.git
```
Then, use
```
cd Spa3D
```
to access the downloaded folder. 

If the "git clone" command does not work with your operation system, you can 
download the zip file from the website 
https://github.com/Lin-Xu-lab/Spa3D.git and decompress it. Then, the folder 
that you need to access is Spa3D-main. 

Note, please see the SPE folder for the processes of spatial pattern enhancement
and 3D data assembly. 

__3.2. Environment setup__

The package has been successuflly tested in a Linux environment of python 
version 3.8.8, pandas version 1.3.4, and g++ version 11.2.0. An option to set up 
the environment is to use Conda 
(https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

You can use the following command to create an environment for Spa3D:
```
conda create -n myenv_Spa3D python=3.8.8 pandas=1.3.4 numpy==1.20.1

pip install numba==0.53.1 python-igraph torch louvain scipy scanpy anndata natsort sklearn
```

When you install the package, if we find a package is not available or is not
consistent with your system, my suggestion is to use the nearest version that 
you could find. For example, if numpy==1.20.1 does not work with your system and
the system messages tells you the earlist version that you could install is 
numpy==1.21.0, please go ahead to install it and see if you could obtain the 
same result with the notebook example below. I have tested numpy==1.21.0 using
my workstation and it works. 

After the environment is created, you can use the following command to activate 
it:
```
conda activate myenv_Spa3D
```

Please install Jupyter Notebook from https://jupyter.org/install. For example, 
you can run
```
pip install notebook
```
in the terminal to install the classic Jupyter Notebook. 

__3.3. Install Spa3D__

To install Spa3D on your computer, please run
```
python setup.py install --user &> log
```
in the terminal.

After doing these successfully, you are supposed to be able to import Spa3D when 
you are using Python or Jupyter Notebook:
```
import Spa3D
```

Another option is to directly use
```
import __init__ as Spa3D
```
to import Spa3D in the current folder. This way is suggested if you would like 
to development your own method based on this Spa3D software.

# 4. Parameter setting

Most of the parameters have a default value. There are two parameters that I 
think are the most important ones:

num_cluster: The number of clusters. 

num_neighbors: The number of nearest neighbors for GCN clustering. 

There are also some other parameters:

p: This parameter adjusts the weights of the current cell and its nearest 
neighbors when searching for the decal value, when p = 0.5, it means their 
weights are equal. 

start: This parameter determines the lowest value of the searching range for the
decay value. 

end: This parameter determines the highest value of the searching range for the 
decay value. 

For introductions about other parameters, please see the *.py files. 


# 5. Examples

__5.1. Spa3D example__

There are one Spa3D examples in the "Spa3D-example" folder.
```
cd Spa3D-example
```
Please use jupyter notebook to open the spagcn_heart3D_week6_20230510.ipynb for 
the 3D human heart example. The input data is a 3D data with spatial pattern 
enhancement. 
  
__5.2. Preprocessing example__

There is one preprocessing example in the "Spa3D-preprocessing-example" folder. 
```
cd preprocessing-example
```
This folder also includes the MPI code for spatial pattern enhancement. 

Please use Jupyter notebook to open the preprocessing_alignment.ipynb for the 
Human Heart example.

__5.3. The notebook script for Spa3D example__
	
Below is the notebook script for the Spa3D example in 5.1. First, please 
decompress the data.zip file and obtain the data folder. Then, please type
```
python3 -m notebook &
```
to open the Jupyter Notebook. Left click the Spa3D_heart_example.ipynb file to 
open it. 

Run the code below to import packages:
```
import sys, os, csv, re, math, cv2, random, torch

import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import numba
import louvain
import scipy
import natsort
import sklearn
import igraph
from scipy.sparse import issparse
from sklearn.metrics.cluster import adjusted_rand_score

import matplotlib.colors as clr
import matplotlib.pyplot as plt
```
Run the code below to import Spa3D:
```
import __init__ as Spa3D
```

Run the code below to print the version of packages:
```
print("Spa3D version: " + Spa3D.__version__)
print("python version: " + sys.version)
print("pandas version: " + pd.__version__)
print("numpy version: " + np.__version__)
print("numba version: " + numba.__version__)
print("python-igraph: " + igraph.__version__)
print("torch version: " + torch.__version__)
print("louvain version: " + louvain.__version__)
print("scipy version: " + scipy.__version__)
print("scanpy version: " + sc.__version__)
print("anndata version: " + ad.__version__)
print("natsort version: " + natsort.__version__)
print("sklearn version: " + sklearn.__version__)
```
The output is:

![Fig](/images/Figure2.png)

Load data and apply simple processings such as norm and log. 
```
adata = sc.read("/home/chentang/data/humanHeart/adata3D_week6_part_aligned_lfft_filter.h5ad")
adata.var_names_make_unique()
Spa3D.preprocessing_filter_by_gene_number(adata,min_cells=3)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
```

Create x_array, y_array, z_array, x_pixel, y_pixel, z_pixel, and set up the
paramters:
```
x_array = adata.obs["x_array"].tolist()
y_array = adata.obs["y_array"].tolist()
z_array = adata.obs["z_array"].tolist()
x_pixel = adata.obs["x_pixel_aligned"].tolist()
y_pixel = adata.obs["y_pixel_aligned"].tolist()
z_pixel = adata.obs["z_pixel"].tolist()

p = 0.5
start = 0.01
end = 200
num_clusters = 5
num_neighbors = 20
r_seed = t_seed = n_seed = 100
```
Calculate the adjacent matrix, decay_value, and the resolution:
```
adj3d = Spa3D.calculate_adj_matrix_3D(x = x_pixel, y = y_pixel, m = z_pixel)
decay_value = Spa3D.search_decay_value_3D(adj3d, p = p, start = start, end = end)
resolution = Spa3D.search_proper_resolution_3D(adata, adj3d, decay_value, num_clusters = num_clusters, num_neighbors = num_neighbors, r_seed = r_seed, t_seed = t_seed, n_seed = n_seed)
```
The output is:

![Fig](/images/Figure3.png)

GCN clustering
```
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)

clf = Spa3D.Spa3D()
clf.set_l(decay_value)
clf.train(adata, adj3d, init_spa = True,init = "louvain",res = resolution, tol=5e-16, lr=0.05, num_neighbors=num_neighbors, max_epochs=2000)
y_pred, prob, pca=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
```
The output is:

![Fig](/images/Figure4.png)

Plot the Spa3D Figure:
```
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]
x_array = adata.obs["x_array"]
y_array = adata.obs["y_array"]
z_array = 12 - adata.obs["z_array"]
z_array_cluster = set(z_array.to_list())
x_pixel = adata.obs["x_pixel_aligned"]
y_pixel = adata.obs["y_pixel_aligned"]
z_pixel = adata.obs["z_pixel"]

fig = plt.figure(figsize=(16, 16))
ax = fig.add_subplot(projection='3d')
ax.grid(False)

dict_x_pixel = {}
dict_y_pixel = {}
dict_z_array = {}
dict_z_pixel = {}
pred = adata.obs["pred"].to_numpy()
cluster = set(pred)
#plt.axis('off')

ax.set_xlim3d(0, 1000)

for j in z_array_cluster:
    
    i = 0

    for id in cluster:
        dict_x_pixel[id] = x_pixel[(pred == id) & (z_array == j)].tolist()
        dict_y_pixel[id] = y_pixel[(pred == id) & (z_array == j)].tolist()
        dict_z_pixel[id] = z_pixel[(pred == id) & (z_array == j)].tolist()
        # ax.axes.invert_yaxis()
        ax.scatter(dict_z_pixel[id], dict_x_pixel[id],dict_y_pixel[id], c=plot_color[i], linewidths = 5.0)
        i = i + 1
```

The output is:
![Fig](/images/Figure5.png)

Plot the Spa3D Figure that can show each individual slide:
```
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]
x_array = adata.obs["x_array"]
y_array = adata.obs["y_array"]
z_array = 12 - adata.obs["z_array"]
z_array_cluster = set(z_array.to_list())
x_pixel = adata.obs["x_pixel_aligned"]
y_pixel = adata.obs["y_pixel_aligned"]
z_pixel = adata.obs["z_pixel"]

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
ax.grid(False)

dict_x_pixel = {}
dict_y_pixel = {}
dict_z_array = {}
pred = adata.obs["pred"].to_numpy()
cluster = set(pred)
plt.axis('off')

for j in z_array_cluster:
    
    i = 0

    for id in cluster:
        dict_x_pixel[id] = x_pixel[(pred == id) & (z_array == j)].tolist()
        dict_y_pixel[id] = y_pixel[(pred == id) & (z_array == j)].tolist()
        dict_z_array[id] = z_array[(pred == id) & (z_array == j)].tolist()
        # ax.axes.invert_yaxis()
        ax.scatter(dict_z_array[id], dict_x_pixel[id],dict_y_pixel[id], c=plot_color[i], linewidths = 2.0)
        i = i + 1
```
The output is:

![Fig](/images/Figure6.png)

Calculate Spa3D ARI:
```
adata.obs["annotation2"] = adata.obs["annotation"]

domains="annotation2"
b = adata.obs[domains].to_numpy()
n = len(b)
for i in range(n):
    if b[i] == 0:
        b[i] = 0
    elif b[i] == 1:
        b[i] = 1
    elif b[i] == 2:
        b[i] = 1
    elif b[i] == 3:
        b[i] = 1
    elif b[i] == 4:
        b[i] = 2
    elif b[i] == 5:
        b[i] = 3
    elif b[i] == 6:
        b[i] = 3
    elif b[i] == 7:
        b[i] = 3
    elif b[i] == 8:
        b[i] = 4
    elif b[i] == 9:
        b[i] = 0
    else:
        b[i] = 0
        
b = b.astype(int)

adata.obs[domains] = b

domains1="pred"
domains2="annotation2"
a = adata.obs[domains1].to_numpy()
b = adata.obs[domains2].to_numpy()
n = len(b)
a = a.astype(float)
b = b.astype(float)
b[np.isnan(b)] = 0

ARI = adjusted_rand_score(a, b)

print("ARI =", ARI)

```
The output is:

![Fig](/images/Figure8.png)

Plot the slide ID 9 from Spa3D:
```
adata_sub = adata[adata.obs['z_array'].isin([9]),:]
#Set colors used
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Spa3D",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

plot_color=["#F56867","#7495D3","#FEB915","#59BE86","#997273"]

#Plot refined spatial domains
domains="annotation2"
adata_sub.obs[domains] = adata_sub.obs[domains].astype(str)
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Reference",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

domains1="pred"
domains2="annotation2"
a = adata_sub.obs[domains1].to_numpy()
b = adata_sub.obs[domains2].to_numpy()
n = len(b)
a = a.astype(float)
b = b.astype(float)
b[np.isnan(b)] = 0

ARI = adjusted_rand_score(a, b)

print("ARI =", ARI)

```
The output is:

![Fig](/images/Figure9.png)

Plot the slide ID 10 from Spa3D::
```
adata_sub = adata[adata.obs['z_array'].isin([10]),:]
#Set colors used
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Spa3D",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

plot_color=["#F56867","#7495D3","#FEB915","#59BE86","#997273"]

#Plot refined spatial domains
domains="annotation2"
adata_sub.obs[domains] = adata_sub.obs[domains].astype(str)
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Reference",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

domains1="pred"
domains2="annotation2"
a = adata_sub.obs[domains1].to_numpy()
b = adata_sub.obs[domains2].to_numpy()
n = len(b)
a = a.astype(float)
b = b.astype(float)
b[np.isnan(b)] = 0

ARI = adjusted_rand_score(a, b)

print("ARI =", ARI)
```
The output is:

![Fig](/images/Figure10.png)

Plot the slide ID 11 from Spa3D::
```
adata_sub = adata[adata.obs['z_array'].isin([11]),:]
#Set colors used
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Spa3D",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

plot_color=["#F56867","#7495D3","#FEB915","#59BE86","#997273"]

#Plot refined spatial domains
domains="annotation2"
adata_sub.obs[domains] = adata_sub.obs[domains].astype(str)
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Reference",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

domains1="pred"
domains2="annotation2"
a = adata_sub.obs[domains1].to_numpy()
b = adata_sub.obs[domains2].to_numpy()
n = len(b)
a = a.astype(float)
b = b.astype(float)
b[np.isnan(b)] = 0

ARI = adjusted_rand_score(a, b)

print("ARI =", ARI)
```
The output is:

![Fig](/images/Figure11.png)

Plot the slide ID 12 from Spa3D::
```
adata_sub = adata[adata.obs['z_array'].isin([12]),:]
#Set colors used
plot_color=["#F56867","#7495D3","#59BE86","#FEB915","#997273"]

#Plot spatial domains
domains="pred"
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Spa3D",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

plot_color=["#F56867","#7495D3","#FEB915","#59BE86","#997273"]

#Plot refined spatial domains
domains="annotation2"
adata_sub.obs[domains] = adata_sub.obs[domains].astype(str)
num_celltype=len(adata_sub.obs[domains].unique())
adata_sub.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata_sub,alpha=1,y="y_pixel_aligned",x="x_pixel_aligned",color=domains,title="Reference",color_map=plot_color,show=False,size=100000/adata_sub.shape[0])
ax.set_aspect('equal', 'box')
#ax.axes.invert_yaxis()
#plt.savefig(pred_pngfilename, dpi=600)
#plt.close()

domains1="pred"
domains2="annotation2"
a = adata_sub.obs[domains1].to_numpy()
b = adata_sub.obs[domains2].to_numpy()
n = len(b)
a = a.astype(float)
b = b.astype(float)
b[np.isnan(b)] = 0

ARI = adjusted_rand_score(a, b)

print("ARI =", ARI)
```
The output is:

![Fig](/images/Figure12.png)

# 6. Contact information

Please contact our team if you have any questions:

Chen Tang (Chen.Tang@UTSouthwestern.edu)

Lei Dong (Lei.Dong@UTSouthwestern.edu)

Xue Xiao (Xiao.Xue@UTSouthwestern.edu)

Yuansheng Zhou (Yuansheng.Zhou@UTSouthwestern.edu)

Lin Xu (Lin.Xu@UTSouthwestern.edu)

Please contact Chen Tang for questions about programming and this README file.

# 7. Copyright information 

The Spa3D software uses the BSD 3-clause license. Please see the "LICENSE" file
for the copyright information. 

Notice: This Spa3D software is adapted from the SpaGCN code 
       (https://github.com/jianhuupenn/SpaGCN). 
       Please see the "LICENSE" file for copyright details of the SpaGCN 
       software. The implementation of the SpaGCN software is described in the 
       publication "SpaGCN: Integrating gene expression, spatial location and 
       histology to identify spatial domains and spatially variable genes by 
       graph convolutional network" 
       (https://doi.org/10.1038/s41592-021-01255-8). 
