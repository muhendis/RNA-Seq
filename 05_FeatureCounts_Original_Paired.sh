#!/bin/bash
# FeatureCount
#PBS -N FeatureCount
#PBS -l nodes=1:regular,walltime=2:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogFC
#PBS -e /scratch/wsspaces/nsharma-test-0/LogFC
#PBS -j oe
#PBS -V

############ Change Directory ############
cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############
module load apps/subread/1.5.0-p3/gcc-4.4.7

###################################### Feature counts ############################
mkdir -p ./FeatureCounts_Original_Paired
cd ./Mapping_STAR
find . -type f -iname "*Aligned.sortedByCoord.out.bam" |
while read filename
do
    fbname=$(basename "$filename" | cut -d. -f1)
    featureCounts -T 16 \
                  -t exon \
                  -s 2 \
                  -p \
                  -g gene_id \
                  -a ../STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94/Homo_sapiens.GRCh38.94.gtf \
                  -o ../FeatureCounts_Original_Paired/${fbname}".fCounts.txt" \
                  $filename 
done

