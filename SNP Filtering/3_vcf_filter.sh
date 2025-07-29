#!/bin/bash
#SBATCH -p eoas_q
#SBATCH -t 4:00:00
#SBATCH -e vcffilt.err
#SBATCH -o vcffilt.log


module load anaconda/3.8.3
eval "$(conda shell.bash hook)"
conda activate vcflib

vcffilter -f "QUAL > 20" COMBINED_n7.vcf  | vcfallelicprimitives >SNP_only_COMBINED_n7.vcf
