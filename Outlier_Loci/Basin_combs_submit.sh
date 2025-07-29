#!/bin/bash

#SBATCH --job-name="Basin_combs"
#SBATCH -e Bsain_combs.err
#SBATCH -N 5
#SBATCH -n 10

module load gnu R/4.1.0

R CMD BATCH --vanilla ./Basin_ind_combs.R Basin_combs.out
