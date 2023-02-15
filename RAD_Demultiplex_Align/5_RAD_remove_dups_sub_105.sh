#!/bin/bash
#Remove PCR duplicates from sorted bam files, then generate flagstats
#run from within directory containing sorted bam files
#SBATCH -p eoas_q
#SBATCH -t 8:00:00
#SBATCH -n 20
#SBATCH -e RAD_rmdups.err
#SBATCH -o RAD_rmdups.log

module load gnu

for file in *_sort.bam

do

sample=`echo $file | cut -f1,2,3,4 -d "_" ` #change to match number of fields in your sample name

echo $sample

samtools rmdup "$sample"_sort.bam ../5_filtered_bam_files/"$sample"_sortfltr.bam

done

cd /gpfs/research/fuenteslab/cc_project/5_filtered_bam_files/
mkdir ./sort_filt_flagstat
for file in *.bam

do

sample=`echo $file | cut -f1 -d "." `

echo $sample
samtools flagstat "$sample".bam \
> ./sort_filt_flagstat/"$sample"_flagstat.txt \
2> ./sort_filt_flagstat/"$sample"_flgstbam.stderr

done


cd ./sort_filt_flagstat/
#this is just combining all the bam stats files for each sample into one text file that is tab delimited so one sample per row so can then easily manipulate to summarize mapping stats etc

for filename in *_sortfltr_flagstat.txt; do
    cat "$filename"|tr '\n' '\t'
    echo "$filename"
done > ./All_RAD_sort_filt_combined.flagstat_tabbed.txt

#few lines to generate a cleaned up summary file of mapping stats post filtering

awk '{print NF}' All_RAD_sort_filt_combined.flagstat_tabbed.txt #88 columns
awk '{print $105 "\t" $1 "\t" $32 "\t" $36 "\t" $47 "\t" $61 "\t" $66 "\t" $69 "\t" $77 "\t" $81 "\t" $84 "\t" $94 }' All_RAD_sort_filt_combined.flagstat_tabbed.txt > All_RAD_sort_filt_combined_flagstat_reformat.txt
sed -i 's/(//g' All_RAD_sort_filt_combined_flagstat_reformat.txt #-i changes it in the original file; could test sending to another outfile first
sed -i 's/%//g' All_RAD_sort_filt_combined_flagstat_reformat.txt
cat ../../scripts_and_keyfiles/short_flagstat_headers All_RAD_sort_filt_combined_flagstat_reformat.txt >All_head_filt_flagstat_reformat.txt
