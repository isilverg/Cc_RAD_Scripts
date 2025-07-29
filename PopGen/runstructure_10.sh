#!/bin/bash

#SBATCH -t 300:00:00
#SBATCH -N 10
#SBATCH -n 60
#SBATCH -e structure_10.err
#SBATCH -o structure_10.log

module load gnu

#mkdir results_f log harvester
#mkdir k1
#mkdir k2
#mkdir k3
#mkdir k4
#mkdir k5
#mkdir k6
#mkdir k7
#mkdir k8
#mkdir k9
#mkdir k10

#cd log
#mkdir k1
#mkdir k2
#mkdir k3
#mkdir k4
#mkdir k5
#mkdir k6
#mkdir k7
#mkdir k8
#mkdir k9
#mkdir k10

#cd ..

cat structureCommands_10 | parallel -j 100%

mv k10  results_f/
#mkdir harvester_input
#cp results_f/k*/*_f harvester_input
#echo 'Your structure run has finished.'
# Run structureHarvester
#./structureHarvester.py --dir harvester_input --out harvester --evanno --clumpp
#echo 'structureHarvester run has finished.'
#Clean up harvester input files.
#zip Combined_cut_Harvester_Upload.zip harvester_input/*
#mv Combined_cut_Harvester_Upload.zip harvester/
#rm -rf harvester_input
