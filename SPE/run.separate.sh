#!/bin/bash
#SBATCH --job-name=gene_denoise_v1
#SBATCH --partition=256GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --output=mutiTaskJob.%j.out
#SBATCH --error=mutiTaskJob.%j.time
#SBATCH --mail-user=chen.tang@utsouthwestern.edu
#SBATCH --mail-type=ALL

module add python/3.8.x-anaconda

python3 -u separate.151673.py
