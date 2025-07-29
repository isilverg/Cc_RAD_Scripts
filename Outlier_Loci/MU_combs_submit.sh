#!/bin/bash

#SBATCH --job-name="MU_combs"
#SBATCH -e MU_combs.err
#SBATCH -p eoas_q
#SBATCH -n 10

module load gnu R/4.1.0

R CMD BATCH --vanilla ./MU_ind_combs.R MU_combs.out
