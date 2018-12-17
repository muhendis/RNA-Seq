#!/bin/bash
# PreFilterQC
#PBS -N PreFilterQC
#PBS -l nodes=1:large,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogPreFilterQC
#PBS -e /scratch/wsspaces/nsharma-test-0/LogPreFilterQC
#PBS -j oe
#PBS -V

############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############

module load apps/fastqc/0.11.3/linux-x86_64
module load apps/multiqc/1.4

############ Create Directories #############

mkdir -p ./QC
mkdir -p ./QC/PreFilteringFastQCzip
mkdir -p ./QC/PreFilteringFastQChtml

############ pre filtering  QC check #############
cd ./Reads 
find -name '*.fastq.gz' | xargs fastqc -t 32  # pass all reads as argument for fastqc
mv *fastqc.zip ../QC/PreFilteringFastQCzip/ 
mv *fastqc.html ../QC/PreFilteringFastQChtml/ 
cd ../QC/PreFilteringFastQCzip
multiqc ./* 
