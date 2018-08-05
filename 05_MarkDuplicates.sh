###################################### Merge lanes ############################
cd ./Mapping_STAR
while read -r line 
do
name=$line
samtools merge -f ./Merged/${name}Aligned.sortedByCoord.out.merged.bam \
               ./${name}_*_L001_R1_001/${name}_*_L001_R1_001Aligned.sortedByCoord.out.bam \
               ./${name}_*_L002_R1_001/${name}_*_L002_R1_001Aligned.sortedByCoord.out.bam 
               #./${name}_*_L003_R1_001/${name}_*_L003_R1_001Aligned.sortedByCoord.out.bam \
               #./${name}_*_L004_R1_001/${name}_*_L004_R1_001Aligned.sortedByCoord.out.bam
done < ../SampleName.txt

################################################### QC_Picardtools ############################
find . -type f -iname '*Aligned.sortedByCoord.out.merged.bam' |
while read filename
do
          fbname=$(basename "$filename" | cut -d. -f1)
          postfix=".sortedByCoord.out.merged.DuplicatesRemoved.bam"
          java -Xmx64g \
          -Djava.io.tmpdir=../tmp \
          -XX:ParallelGCThreads=26 \
          -jar /home/nsharma/anaconda3/share/picard-2.18.9-0/picard.jar \
          MarkDuplicates \
          REMOVE_DUPLICATES=true \
          INPUT=$filename \
          OUTPUT=./Duplicates_Removed/${fbname}${postfix} \
          METRICS_FILE=./Duplicates_Removed/${fbname}_mark_dups_metrics.txt
done
rm -rf ../tmp
cd ../