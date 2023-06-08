import sys
import pandas as pd
import numpy as np

def get_coords(index):
    coords = pd.DataFrame(index=index)
    coords['x'] = index.str.split('x').str.get(0).map(float)
    coords['y'] = index.str.split('x').str.get(1).map(float)
    return coords
    
def main():

	counts = pd.read_csv('../data/adata_heart_9_data.csv', index_col = 0)
	geneList = counts.columns.to_numpy()
	np.savetxt("../data/geneList_9.out", geneList, delimiter=',', fmt = '%s')

	i = 0

	for gene in geneList:
		sample_info = get_coords(counts.index)
		sample_info[gene] = counts[gene]
		filename = gene + ".csv"
		folderpath = "../data/heart9/"
		file_location = folderpath + filename
		print(i, file_location)
		sample_info.to_csv(file_location)
			
		i = i + 1
		
if __name__ == "__main__":
	sys.exit(main())	
