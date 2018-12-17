#!/bin/bash
# NGS
#PBS -N NGS
#PBS -l nodes=1:large,walltime=12:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogNGS
#PBS -e /scratch/wsspaces/nsharma-test-0/LogNGS
#PBS -j oe
#PBS -V


############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############

module load apps/fastqc/0.11.3/linux-x86_64
module load apps/multiqc/1.4
module load apps/bbmap/36.20
module load apps/star/2.5.1b/gcc-5.1.0
module load apps/samtools/1.3.1/gcc-4.4.7
module load apps/subread/1.5.0-p3/gcc-4.4.7


############ Create Directories #############

mkdir -p ./QC
mkdir -p ./QC/PreFilteringFastQCzip
mkdir -p ./QC/PreFilteringFastQChtml
mkdir -p ./QC/PostFilteringFastQCzip
mkdir -p ./QC/PostFilteringFastQChtml
mkdir -p ./QC/QCpassed
mkdir -p ./QC/QCfailed
mkdir -p ./Reads/temp_unmerged
mkdir -p ./STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94
mkdir -p ./Mapping_STAR
mkdir -p ./Mapping_STAR/multiqc
mkdir -p ./tmp
mkdir -p ./FeatureCounts_Original_Paired

############## Main Script #######

qsub 01_PreFilteringQC.sh
qsub 02_Filtering.sh
qsub 03_StarCreateGenomeIndices.sh
qsub 04_Mapping.sh
qsub 05_FeatureCounts_Original_Paired.sh
R CMD 06_DESEQ2.R
