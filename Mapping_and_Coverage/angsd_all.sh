#!/bin/bash -l

#SBATCH -p eoas_q
#SBATCH -N 4
#SBATCH -t 40:00:00
#SBATCH -o angsd_60.out
#SBATCH -e angsd_60.err
#SBATCH --job-name="Cc_60"

nInd=94

minInd=60

module load gnu openmpi
angsd -bam bamlist_all -out all_60 -minQ 20 -minMapQ 10 -minInd $minInd -GL 1 -doMajorMinor 1 -doMaf 2
