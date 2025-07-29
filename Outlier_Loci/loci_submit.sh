#!/bin/bash

#SBATCH --job-name="Loc_MU"
#SBATCH -e Loc_MU.err
#SBATCH --mem=80000M

module load gnu R/4.1.0

R CMD BATCH --vanilla ./loci_process_10_2_24.R Loc_MU.out

