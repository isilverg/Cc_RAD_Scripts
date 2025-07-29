#!/bin/bash

#SBATCH -n 4
#SBATCH -p eoas_q
#SBATCH -t 4:00:00
#SBATCH -e admixture.err

pops=10

for K in $( seq 0 $pops )

do

./admixture --cv cc_combined.bed $K | tee log${K}.out

done
