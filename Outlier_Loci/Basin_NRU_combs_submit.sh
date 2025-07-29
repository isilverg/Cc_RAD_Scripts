#!/bin/bash

#SBATCH --job-name="Basin_NRU_combs"
#SBATCH -e Bsain_combs.err
#SBATCH -N 1
#SBATCH -n 10

module load gnu R/4.1.0

R CMD BATCH --vanilla ./Basin_NRU_ind_combs.R Basin_NRU_combs.out
