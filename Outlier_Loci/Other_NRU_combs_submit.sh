#!/bin/bash

#SBATCH --job-name="Other_NRU_combs"
#SBATCH -e Other_NRU_combs.err
#SBATCH -N 1
#SBATCH -n 10

module load gnu R/4.1.0

R CMD BATCH --vanilla ./Other_NRU_ind_combs.R Other_NRU_combs.out
