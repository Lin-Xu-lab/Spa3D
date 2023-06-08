#!/bin/bash
#SBATCH --job-name=gene_combine
#SBATCH --partition=256GBv1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --output=mutiTaskJob.%j.out
#SBATCH --error=mutiTaskJob.%j.time
#SBATCH --mail-user=chen.tang@utsouthwestern.edu
#SBATCH --mail-type=ALL

module add python/3.8.x-anaconda

python3 -u combine.py
