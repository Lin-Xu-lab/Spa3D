from mpi4py import MPI
import sys
import __main__
import math 
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import pywt
from scipy.signal import hilbert, hilbert2

def get_coords(index):

    coords = pd.DataFrame(index=index)
    coords['x'] = index.str.split('x').str.get(0).map(float)
    coords['y'] = index.str.split('x').str.get(1).map(float)
    
    return coords
 
def distribute(p, src_begin, src_end, start_point, length):
    n = src_end - src_begin + 1
    len = int(n / p)
    nre = n % p
    len_max = len
    
    for i in range(p):
        length[i] = len
        
    if nre != 0:
        for i in range(nre):
            length[i] += 1
        len_max += 1
        
    start_point[0] = src_begin
    
    for i in range(1, p + 1):
        start_point[i] = start_point[i - 1] + length[i - 1]
    
    return len_max  
    
def spe_alft(rank, gene, sample_info, sample_info_spe_alft):

	sample_info_position = get_coords(sample_info.index)
	n = len(sample_info_position)
		
	para_neighbors = 25
	para_neighbors_threashold = 17
	dx = 30.0
	dy = 30.0

	nbrs = NearestNeighbors(n_neighbors = para_neighbors, algorithm='ball_tree').fit(sample_info_position)
	distances, indices = nbrs.kneighbors(sample_info_position)

	for i in range(n):
		sample_info_local = sample_info.iloc[indices[i]]
		
		x_min = int(sample_info_local['x'].min())
		y_min = int(sample_info_local['y'].min())
		x_max = int(sample_info_local['x'].max())
		y_max = int(sample_info_local['y'].max())
		x_range = int((x_max - x_min) / dx + 0.5) + 2
		y_range = int((y_max - y_min) / dy + 0.5) + 2
		
		a = np.zeros((y_range, x_range) , dtype = np.float32)
		
		for j in range(para_neighbors_threashold):   
			x = int((sample_info['x'][indices[i][j]] - x_min) / dx + 0.5)
			y = int((sample_info['y'][indices[i][j]] - y_min) / dy + 0.5)
			a[y][x] = sample_info[gene][indices[i][j]]
#-------------------------
		a_fft2 = np.fft.rfft2(a)
		a_fft2_abs = abs(a_fft2) 	
		sum_energy = a_fft2_abs.sum()
		
		ny = np.shape(a_fft2)[0]
		nx = np.shape(a_fft2)[1]	        
		a_fft2_error = np.zeros((ny, nx) , dtype = np.complex128)
			
		perc = 0.02
		sum_energy_iter = 0.0
		a_fft2_abs_index_rank = np.dstack(np.unravel_index(np.argsort(a_fft2_abs.ravel()), a_fft2_abs.shape))[0]
		
		if sum_energy * perc > 1.0e-16:
			iter_alft = 0		  
			while sum_energy_iter < sum_energy * perc:
				max_index = a_fft2_abs_index_rank[iter_alft]
				iy_max = max_index[0]
				ix_max = max_index[1]

				a_fft2_error[iy_max][ix_max] = a_fft2[iy_max][ix_max]
				sum_energy_iter += a_fft2_abs[iy_max][ix_max]
				
				iter_alft += 1

		a_fft2_reserved = a_fft2 - a_fft2_error

		a_inv = np.fft.irfft2(a_fft2_reserved)
#-------------------------
		x0 = int((sample_info['x'][indices[i][0]] - x_min) / dx + 0.5)
		y0 = int((sample_info['y'][indices[i][0]] - y_min) / dy + 0.5)
		sample_info_spe_alft[gene][indices[i][0]] = a_inv[y0][x0]

	if rank == 0:
		print ("para_neighbors = ", para_neighbors, "para_neighbors_threashold = ", para_neighbors_threashold, "dx = ", dx, "dy = ", dy, "perc = ", perc)
    
def main():
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	p = comm.Get_size()
	
	if rank == 0:
		print ("p =", p)
	sys.stdout.flush()

	for id in range(9, 13):
		filename_genelist = "data/geneList_" + str(id) + ".out"
		
		if rank == 0:
			geneList = np.loadtxt(filename_genelist, delimiter=',', dtype = str)
		else:
			geneList = None
			
		geneList = comm.bcast(geneList, root = 0)
		
		src_begin = 0
		src_end = len(geneList) - 1
		
		if rank == 0:
			print ("src_begin = ", src_begin, "src_end = ", src_end)
			sys.stdout.flush()

		start_point = np.zeros((p + 1, ), dtype = int)
		length = np.zeros((p, ), dtype = int)
		
		len_max = distribute(p, src_begin, src_end, start_point, length)

		if rank % 100 == 0:	
			print ("rank = ", rank, "start = ", start_point[rank], "end = ", start_point[rank + 1] - 1, "len = ", length[rank])
			sys.stdout.flush()
		
		ig = 0
		
		for isrc in range(start_point[rank], start_point[rank + 1]):
			gene = geneList[isrc]
			
			filename = gene + ".csv"
			folderpath = "data/heart" + str(id) + "/"
			file_location = folderpath + filename
			if rank % 100 == 0:
				print("rank = ", rank, "input = ", file_location)
				sys.stdout.flush()
				
			sample_info = pd.read_csv(file_location, index_col = 0)

			if rank == 0:
				print ("rank = ", rank, "ig = ", ig, "gene = ", gene)
				sys.stdout.flush()
				ig = ig + 1

			sample_info_spe_alft = sample_info.copy()
			spe_alft(rank, gene, sample_info, sample_info_spe_alft)
			
			filename = gene + ".csv"
			filepath = "data/heart_out_" + str(id) + "/"
			file_location = filepath + filename
			if rank % 100 == 0:
				print("rank = ", rank, "out = ", file_location)
				sys.stdout.flush()
			sample_info_spe_alft.to_csv(file_location)
		
if __name__ == "__main__":
	sys.exit(main())	
