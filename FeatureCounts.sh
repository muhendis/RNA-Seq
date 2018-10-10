#!/bin/bash
# FCount
#PBS -N FCount
#PBS -l nodes=1:ppn=16,mem=32Gb,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCountsDeduped
#PBS -e /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCountsDeduped
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Change Directory ############

cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/Mapping_STAR/Duplicates_Removed

###################################### Feature counts ############################
mkdir -p ../FeatureCounts
find . -type f -iname '*Aligned.sortedByCoord.out.Deduped.bam' |
while read filename
do
    fbname=$(basename "$filename" | cut -d. -f1)
    featureCounts -T 16 \
                  -t exon \
                  -s 1 \
                  -g gene_id \
                  -a ../../GenomeIndicesSTAR_hg38/hg38.gtf \
                  -o ../FeatureCounts_Duplicates_Removed/${fbname}".fCounts.txt" \
                  $filename 
done