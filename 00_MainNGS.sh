###########################################################################################################################
################################     	   Script to analyse RNA-Seq Data               ###################################
################################        Created by Nitin Sharma 2  July 2018            ###################################
###########################################################################################################################

#!/bin/bash
# NGS
#PBS -N NGS
#PBS -l nodes=1:ppn=32,mem=2000Gb,walltime=10:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/RNA-Seq_Mouse/LogNGS
#PBS -e /scratch/wsspaces/nsharma-test-0/RNA-Seq_Mouse/LogNGS
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/RNA-Seq_Mouse/
  
  ############ Load modules ############

module load apps/bbmap/36.20
module load apps/fastqc/0.11.3/linux-x86_64
module load apps/multiqc/1.4
module load apps/star/2.5.1b/gcc-5.1.0
module load apps/picardtools/1.96/noarch
module load apps/samtools/1.3.1/gcc-4.4.7
module load apps/subread/1.5.0-p3/gcc-4.4.7

############ Create Directories #############

mkdir -p ./QC
mkdir -p ./QC/PreFilteringFastQCzip_merged
mkdir -p ./QC/PreFilteringFastQChtml_merged
mkdir -p ./QC/PostFilteringFastQCzip_merged
mkdir -p ./QC/PostFilteringFastQChtml_merged
mkdir -p ./QC/QCpassed
mkdir -p ./QC/QCfailed
mkdir -p ./Reads/temp_unmerged
mkdir -p ./Reads/temp_merged_prefiltering
mkdir -p ./Reads/temp_merged_postfiltering
mkdir -p ./Mapping_STAR
mkdir -p ./Mapping_STAR/Merged
mkdir -p ./Mapping_STAR/Duplicates_Removed
mkdir -p ./tmp
mkdir -p ./MarkedDuplicates
mkdir -p ./GenomeIndicesGRCm38.93_STAR
mkdir -p ./FeatureCounts

############## Source files #######

./01_PreFilteringQC.sh
./02_FilteringQC.sh
./03_PostFilteringQC.sh
./04_Mapping.sh
./05_MarkDuplicates.sh
./06_FeatureCounts.sh
