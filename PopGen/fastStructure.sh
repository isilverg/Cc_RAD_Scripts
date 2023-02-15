#!/bin/bash

#SBATCH -p eoas_q
#SBATCH -t 1:00:20
#SBATCH -n 20
#SBATCH -e fastStructure.err
#SBATCH -o fastStructure.log

inputfile=~/vcftools_out/CCGOM_raw_filtered.vcf
fileprefix=CCGOM_filter
outfileprefix=~/CCGOM_fastStructure/CCGOM

# convert vcf to plink bed
#~/plink/plink --vcf $inputfile --double-id --allow-extra-chr --out ~/CCGOM_fastStructure/$fileprefix

#activate environment
source ~/fastStructure/bin/activate

# run fastStructure
python ~/proj/fastStructure/structure.py -K 1 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 2 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 3 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 4 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 5 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 6 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 7 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 8 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 9 --input=./$fileprefix --output=$outfileprefix
python ~/proj/fastStructure/structure.py -K 10 --input=./$fileprefix --output=$outfileprefix

#run chooseK
python ~/proj/fastStructure/chooseK.py --input=$outfileprefix > ${outfileprefix}chooseKout.txt

#extract marginal likelihood
touch ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.1.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.2.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.3.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.4.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.5.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.6.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.7.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.8.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.9.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
grep 'Marginal Likelihood =' ~/CCGOM_fastStructure/CCGOM.10.log | cut -d " " -f 4 >> ${outfileprefix}_ML.txt
