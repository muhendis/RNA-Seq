#!/bin/bash
# FCount
#PBS -N FCount
#PBS -l nodes=1:ppn=16,mem=32Gb,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCounts
#PBS -e /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCounts
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Load modules ############
module load apps/multiqc/1.4

############ Change Directory ############

cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/
  
############ run Multiqc ############
cd ./Reads/temp_merged_prefiltering
mv *fastqc.zip ../../QC/PreFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv *fastqc.html ../../QC/PreFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
cd ../../QC/PreFilteringFastQCzip_merged 
multiqc ./* 

cd ../../Reads/temp_merged_postfiltering # change directory to temp_merged_prefiltering
mv *fastqc.zip ../../QC/PostFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv *fastqc.html ../../QC/PostFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
cd ../../QC/PostFilteringFastQCzip_merged 
multiqc ./* 

cd ../../Mapping_STAR
mkdir -p ./multiqc
find . -type f -iname "[a-zA-Z0-9]*Log.final.out" | xargs cp -t ./multiqc
find . -type f -iname "[a-zA-Z0-9]*ReadsPerGene.out.tab" | xargs cp -t ./multiqc
cd ./multiqc
multiqc ./* 
