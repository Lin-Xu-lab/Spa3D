#!/bin/bash
#SBATCH --job-name=gene_denoise_v1
#SBATCH --partition=256GBv1
#SBATCH --nodes=8
#SBATCH --ntasks=200
#SBATCH --time=24:00:00
#SBATCH --output=mutiTaskJob.%j.out
#SBATCH --error=mutiTaskJob.%j.time
#SBATCH --mail-user=chen.tang@utsouthwestern.edu
#SBATCH --mail-type=ALL

module add python/3.8.x-anaconda
module add mpich/ge/gcc/64/3.3.2
module add mpich/ge/gcc/64/3.2rc2
module add mpiexec/0.84_432

mpirun python3 mpi.151673.lfft.new.py
