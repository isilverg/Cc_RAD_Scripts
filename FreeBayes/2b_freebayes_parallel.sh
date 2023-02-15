#!/bin/bash

#use freebayes to generate a vcf file with all samples in a dataset

#SBATCH -t 60:00:00
#SBATCH -n 60
#SBATCH -p eoas_q
#SBATCH -e freebayes-p.err
#SBATCH -o freebayes-p.log

module load gnu
module load anaconda/3.8.3
eval "$(conda shell.bash hook)"
conda activate /gpfs/home/isilvergorges/.conda/envs/FreeBayes/

reference=../../reference_genome/GCF_02*.fna

freebayes-parallel ./targets.bed 60 -f $reference -L ./bamlist.txt --use-best-n-alleles 7 > CCGOM_n7.vcf

conda deactivate && conda deactivate
