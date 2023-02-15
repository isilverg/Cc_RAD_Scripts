#!/bin/bash

#use samtools index to prepare files for freebayes to generate a vcf file with all samples in a dataset

#SBATCH -t 60:00:00
#SBATCH -n 60
#SBATCH -p eoas_q
#SBATCH -e samindex.err
#SBATCH -o samindex.log

module load gnu

#first add readgroup ids to bam files in case they aren't there already with picard- if unsure, just do this step
#run this script from within the directory containing your sorted, filtered bam files

reference=../../reference_genome/GCF_02*.fna

for file in ../5*/*sortfltr.bam
do
date
echo $file
sample=$(basename $file _sortfltr.bam )
#echo $sample
java -jar ~/software/picard/picard.jar AddOrReplaceReadGroups  I= $file O= ./"$sample"_RG_sortfltr.bam RGID= $sample RGLB= $sample RGPL=illumina RGPU=unit1 RGSM= $sample
done

#index new bam files
for file in ./*_RG_sortfltr.bam
do
samtools index $file
done

#make new bam file list
ls ./*_RG_sortfltr.bam > bamlist.txt

