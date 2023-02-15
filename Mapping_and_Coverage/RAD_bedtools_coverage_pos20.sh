#!/bin/bash
# run bedtools coverage for file in directory *.bam

#SBATCH -p eoas_q
#SBATCH -t 2:00:00
#SBATCH -n 15
#SBATCH -e bedtools_pos20.err
#SBATCH -o bedtools_pos20.log

#module load bedtools/2.28.0

#mkdir ../LMK_coverage_scratch/bedtools_pos20_cov

for file in ../5*/*.bam

do

sample=`echo $file | cut -d "_" -f 4-7 | cut -d "/" -f 2 ` \

echo $sample
~/bedtools/bedtools coverage -a ./RAD_bedtools_pos20_134_60.bed -b ../5*/"$sample"_sortfltr.bam  >../cc_test_coverage/bedtools_pos20_60_cov/"$sample".pos20.cov.txt

done
