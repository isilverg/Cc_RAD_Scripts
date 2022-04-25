#!/bin/bash -l

#SBATCH -p eoas_q
#SBATCH -N 10
#SBATCH -t 8:00:00
#SBATCH -o angsd.out
#SBATCH -e angsd.err

nInd=96

minInd=1

module load gnu openmpi
angsd -bam bamlist_all -out all -minQ 20 -minMapQ 10 -minInd $minInd -GL 1 -doMajorMinor 1 -doMaf 2 
