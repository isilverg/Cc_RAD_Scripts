#!/bin/bash
#SBATCH -p eoas_q
#SBATCH -t 40:00:00
#SBATCH -n 20
#SBATCH -e depth.err
#SBATCH -o depth.log

#generate a region list for freebayes from one Bam file breaking apart by depth into 10000 regions

module load gnu
module load anaconda/3.8.3
eval "$(conda shell.bash hook)"
conda activate FreeBayes

samtools faidx ../../reference_genome/GCF_02*.fna #make sure to index your ref genome to produce an .fna.fai file
samtools depth ../5*/*wG04*_sortfltr.bam | ./coverage_to_regions.py ../../reference_genome/GCF_02*.fna.fai 10000 >targets.bed

