#!/bin/bash
# FCount
#PBS -N FCount
#PBS -l nodes=1:ppn=16,mem=32Gb,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCountsKnownGeneDedupedPaired
#PBS -e /scratch/wsspaces/nsharma-Alessio_ACseries-0/LogFCountsKnownGeneDedupedPaired
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Load modules ############
module load apps/subread/1.5.0-p3/gcc-4.4.7

############ Change Directory ############
cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/Mapping_STAR/Duplicates_Removed

###################################### Feature counts ############################
mkdir -p ../FeatureCounts_Deduped_Paired
find . -type f -iname "*Aligned.sortedByCoord.out.Deduped.bam" |
while read filename
do
    fbname=$(basename "$filename" | cut -d. -f1)
    featureCounts -T 16 \
                  -t exon \
                  -s 1 \
                  -p \
                  -g gene_id \
                  -a ../../GenomeIndicesSTAR_hg38/hg38.gtf \
                  -o ../FeatureCounts_Deduped_Paired/${fbname}".fCounts.txt" \
                  $filename 
done