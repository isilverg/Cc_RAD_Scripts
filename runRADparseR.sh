#!/bin/bash -l

#SBATCH -t 40:00:00
#SBATCH -n 20
#SBATCH -e all_50.err
#SBATCH -o all_50.log
#SBATCH -p eoas_q

module load R

R CMD BATCH my_RAD_parse_test_exact_length.R
