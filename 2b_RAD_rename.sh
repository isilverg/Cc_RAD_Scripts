#!/bin/bash -l

list=$1
#enter metadata file as argument 1 in command line (cc_keyfile.txt)
x=1
while [ $x -le 96 ]
do

	string="sed -n ${x}p ${list}"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2, $3, $4}')
	set -- $var
	c1=$1 #well barcode
	c2=$2 #RAD well ID
	c3=$3 #SampleID
	c4=$4 #Location

#run with just RA and RB files in different subfolders; run for RA, then edit accordingly and run for RB
#run with echo to check is correct, then if all looks good remove and run without.
mv ../1_raw_fastq/CC_RB/RAD-CCGOM_RB_GG${c1}TGCAGG.fastq ../2_demultiplexed_renamed_fastq/RAD-CCGOM_${c2}_${c3}_${c4}_RB.fastq
#RAD-CCGOM_RA_GGAACCGAGATGCAGG.fastq wA01_MML257_NGM_RA.fastq
	x=$(( $x + 1 ))

done
