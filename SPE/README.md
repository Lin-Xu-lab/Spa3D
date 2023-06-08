The spatial pattern enhancement (SPE) and 3D data assembly processes involve a 
few steps. You also need to decompress the data.zip for your test. 

First, you need to use separate.py and run.slurm.separate.sh to separate the 
slide into different files at different genes. This needs to be done to all the 
four slides with an ID from 9 to 12. 

Second, you need to use mpi.spe.alft.py and run.slurm.spe.sh to do the SPE for 
all the gene folders from 9 to 12. mpi4py needs to be installed for this step. 

Third, you need to use combine.py and run.combine.slurm.sh to combine the gene
folder into SPE slides with IDs from 9 to 12. 

Fourth, you need to use test-heart-create2Danndata.ipynb to convert the SPE csv 
files into h5ad files. 

Fifth, you need to use test-heart-create3Danndata-alignment.ipynb to do the 
pairwise alignment and 3D data assembly. During this process, you need to 
install https://github.com/raphael-group/paste for the pairwise alignment. 
