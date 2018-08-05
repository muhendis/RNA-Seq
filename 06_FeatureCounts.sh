###################################### Feature counts ############################
cd ./Mapping_STAR/Duplicates_Removed
find . -type f -iname '*Aligned.sortedByCoord.out.merged.DuplicatesRemoved.bam' |
while read filename
do
	#fbname=$(basename "$filename" | cut -d. -f1)
	featureCounts \
				-T 32 \
				-t exon \
				-s 1 \
				-g gene_id \
	 			-a ../../GenomeIndicesGRCm38.93_STAR/Mus_musculus.GRCm38.93.gtf \
	 			-o ../../FeatureCounts/${filename}.txt ${filename}
done