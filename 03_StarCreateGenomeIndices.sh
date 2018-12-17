#!/bin/bash
# STARindex
#PBS -N STARindex
#PBS -l nodes=1:large,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogSTAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh3894
#PBS -e /scratch/wsspaces/nsharma-test-0/LogSTAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh3894
#PBS -j oe
#PBS -V

############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############

module load apps/star/2.5.1b/gcc-5.1.0

############ Create Directories #############
mkdir -p ./STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94

############ Run Command ############

indexDir="/scratch/wsspaces/nsharma-test-0/STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94"

STAR --runThreadN 25  \
     --runMode genomeGenerate \
     --genomeDir $indexDir \
     --genomeFastaFiles $indexDir/Homo_sapiens.GRCh38.94.fa \
     --sjdbGTFfile $indexDir/Homo_sapiens.GRCh38.94.gtf \
     --sjdbOverhang 100