#!/bin/bash

#SBATCH -n 4
#SBATCH -p eoas_q
#SBATCH -t 1:00:00
#SBATCH --mail-type=ALL
#SBATCH -o multiqc.out

module load python
import multiqc
multiqc ./ --interactive
