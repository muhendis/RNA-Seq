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

############ Change Directory ############

cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/
  
###################################### Feature counts ############################

############ run Multiqc ############
cd ./QC/PreFilteringFastQCzip_merged 
multiqc ./* 
cd ../PostFilteringFastQCzip_merged 
multiqc ./* 
cd ../../Mapping_STAR/multiqc
multiqc ./* 
rm -rf ../tmp
#rm -rf ./Reads/temp_merged_prefiltering # delete folder temp_merged_prefiltering as that was only for QC check
#rm -rf ./Reads/temp_unmerged # delete folder temp_unmerged as that was only for QC check
#rm -rf ./Reads/temp_merged_postfiltering # delete folder temp_merged_postfiltering as that was only for QC check