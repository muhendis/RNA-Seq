#!/bin/bash
# FeatureCountDeduped
#PBS -N FeatureCountDeduped
#PBS -l nodes=1:regular,walltime=2:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogFCDeduped
#PBS -e /scratch/wsspaces/nsharma-test-0/LogFCDeduped
#PBS -j oe
#PBS -V

############ Change Directory ############
cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############
module load apps/subread/1.5.0-p3/gcc-4.4.7

###################################### Feature counts ############################
mkdir -p ./FeatureCounts_Deduped_Paired
cd ./Mapping_STAR/Duplicates_Removed
find . -type f -iname "*Aligned.sortedByCoord.out.Deduped.bam" |
while read filename
do
    fbname=$(basename "$filename" | cut -d. -f1)
    featureCounts -T 16 \
                  -t exon \
                  -s 1 \
                  -p \
                  -g gene_id \
                  -a ../STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94/Homo_sapiens.GRCh38.94.gtf \
                  -o ../FeatureCounts_Deduped_Paired/${fbname}".fCounts.txt" \
                  $filename 
done