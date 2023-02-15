#!/bin/bash
#run fastqc

#SBATCH -n 4
#SBATCH -p eoas_q
#SBATCH -t 4:00:00
#SBATCH -o fastqcRA.out

module load fastqc

for file in ./*RB.fastq

do

echo $file

fastqc $file

done

