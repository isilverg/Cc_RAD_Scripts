#!/bin/bash

#SBATCH --job-name="MUs_2_gen"
#SBATCH -e GSI_MUs_2-502.err
#SBATCH -n 56
#SBATCH --mem-per-cpu=8000M

module load gnu R/4.1.0

R CMD BATCH --vanilla ./GSI_Cluster_Parallel_MUs_6_3_24.R GSI_MU_2-502.out
