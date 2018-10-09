#!/bin/bash
# GenomeIndicesSTARhg38
#PBS -N GenomeIndicesSTARhg38
#PBS -l nodes=1:ppn=32,mem=90154756kb,walltime=4:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogGenomeIndicesSTARhg38
#PBS -e /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogGenomeIndicesSTARhg38
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Change Directory ############

cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/

############ Load modules ############

module load apps/star/2.5.1b/gcc-5.1.0

############ Run Command ############

STAR --runThreadN 32  --runMode genomeGenerate \
     --genomeDir /scratch/wsspaces/nsharma-Alessio_ACseries-0/GenomeIndicesSTAR_hg38_b \
     --genomeFastaFiles /scratch/wsspaces/nsharma-Alessio_ACseries-0/GenomeIndicesSTAR_hg38_b/hg38.fa \
     --sjdbGTFfile /scratch/wsspaces/nsharma-Alessio_ACseries-0/GenomeIndicesSTAR_hg38_b/hg38.gtf