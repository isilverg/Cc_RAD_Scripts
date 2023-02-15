#!/bin/bash

#SBATCH -o filter_VCF_SNP.log
#SBATCH -e filter_VCF_SNP.err
#SBATCH -n 8
#SBATCH -t 05:00:00
#SBATCH -p eoas_q

module load gnu

VCFTOOLSIMG=~/vcftools.sif
input_vcf=./SNP_only_CCGOM_n7.vcf
fileprefix=CCGOM_raw

##############################################
### Filter VCF for use in pop gen analyses ###
##############################################

#--min-meanDP and --max-meanDP filter by the mean minimum and maximum depth values across all included individuals.
#--maf is minor allele frequency
#--hwe - check for hwe and remove sites with a p value below this threshold
#--max-missing removes sites on the proportion of missing data, where 0 means sites that are totally missing are fine, and 1 means no missing data allowed.
#--recode outputs a new vcf file after filtering

#First filtering step
singularity exec $VCFTOOLSIMG vcftools --vcf $input_vcf --max-missing 0.5 --min-meanDP 10 --max-meanDP 1000 --recode --mac 3 --stdout | bcftools +prune -l 0.5 -w 100000 - -Ov -o ./vcftools_out/CCGOM_raw.vcf
#Second step calculates data missingness
## Per individual ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --missing-indv

## Per site ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --missing-site

#Fourth step calculates deviation from HWE
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --hardy

#Fifth step calculates depth per individual and locus
## Mean depth per individual ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --depth

## mean depth per site averaged across individuals ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --site-mean-depth

## Depth for each genotype in the VCF file ##
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --geno-depth

# Sixth step calculates heterozygosity per individual
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --het

# Seventh step grabs SNP call quality
singularity exec $VCFTOOLSIMG vcftools --vcf ./vcftools_out/${fileprefix}.vcf --out ./vcftools_out/${fileprefix} --site-quality
