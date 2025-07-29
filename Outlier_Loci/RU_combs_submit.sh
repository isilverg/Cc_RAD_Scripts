#!/bin/bash

#SBATCH --job-name="RU_combs"
#SBATCH -e RU_combs.err
#SBATCH -N 10
#SBATCH -n 10

module load gnu R/4.1.0

R CMD BATCH --vanilla ./RU_ind_combs.R RU_combs.out
