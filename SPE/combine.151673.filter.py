import sys
import pandas as pd
import numpy as np

def get_coords(index):
    coords = pd.DataFrame(index=index)
    coords['x'] = index.str.split('x').str.get(0).map(float)
    coords['y'] = index.str.split('x').str.get(1).map(float)
    return coords
    
def main():

	for id in range(17, 18):

		filename_data_input = "/work/PCDC/s182862/data/heart/adata_heart_" + str(id) + "_data.csv"
		filename_data_output = "/work/PCDC/s182862/data/heart/adata_heart_" + str(id) + "_data_lfft_filter_25_17_30x30_002.csv"

		counts = pd.read_csv(filename_data_input, index_col = 0)
		geneList = counts.columns.to_numpy()

		counts_denoise = counts.copy()
		i = 0
		for gene in geneList:
			filename = gene + ".csv"
			folderpath = "/work/PCDC/s182862/data/heart/heart_out_lfft_" + str(id) + "/"
			file_location = folderpath + filename
			print (" i = ", i, "file = ", file_location)
			sample_info = pd.read_csv(file_location, index_col = 0)
			counts_denoise[gene] = sample_info[gene]
			i = i + 1
			
		counts_denoise.to_csv(filename_data_output)
		counts_max = counts_denoise.max().max()
		counts_min = counts_denoise.min().min()
		print(counts_min, counts_max)
		
if __name__ == "__main__":
	sys.exit(main())	
