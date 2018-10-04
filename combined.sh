# reads sample names from SampleName.txt 
cd ./Reads 
while read -r line 
do
name="$line"
############ Merge lanes of the samples for QC check  ############
cat ${name}*.fastq.gz > ./temp_merged_prefiltering/${name}.combined.fastq.gz
############ merged pre filtering  QC check #############
cd ./temp_merged_prefiltering 
find -name ${name}.combined.fastq.gz | xargs fastqc -t 16 
mv ${name}*fastqc.zip ../../QC/PreFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv ${name}*fastqc.html ../../QC/PreFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
# unzip file to "temp_unmerged"
cd ../Reads
find . -type f -iname '${name}*.fastq.gz' | while read filename
do
fbname=$(basename "$filename")
s1=`basename "$filename" | cut -d . -f -2`
gunzip -c -N ${filename}  > ./temp_unmerged/$s1
done
############ run BBDUk for trimming/Filtering on unmerged files ############
cd ./temp_unmerged  # change directory to "temp_unmerged" to run QC check using BBDUK 
find . -type f -iname '${name}*.fastq' | while read filename
do
	fbname=$(basename "$filename")
 	bbduk.sh in=$filename \
  out=../../QC/QCpassed/${fbname}.gz \
	outm=../../QC/QCfailed/${fbname}.gz \
	ref=/apps/modules/pkg/apps/bbmap/36.20/resources/adapters.fa \
	overwrite=t \
	forcetrimleft=11 \
	ktrim=r k=13 \
	mink=5 \
	qtrim=rl trimq=10 \
	hdist=1 \
	minlength=20 \
	stats=../../QC/bbduk${fbaname}.stats.txt
done 
cd ../../ 
############ merged  QC passed reaads #############
cd ./QC/QCpassed # change directory to QCpassed
# reads sample names from SampleName.txt and merge the files with output in folder "temp_merged_prefiltering"
cat ${name}*.fastq.gz > ../../Reads/temp_merged_postfiltering/${name}.QCpassed.combined.fastq.gz
########### FastQC for merged QC passed reads ############
cd ../../Reads/temp_merged_postfiltering # change directory to temp_merged_prefiltering
find -name '${name}.QCpassed.combined.fastq.gz' | xargs fastqc -t 16 # pass all merged reads as argument for fastqc
mv ${name}*fastqc.zip ../../QC/PostFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv ${name}*fastqc.html ../../QC/PostFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
cd ../../

############ Mapping By STAR ############
cd ./QC/QCpassed/
find . -type f -iname '${name}*.fastq.gz' |
while read filename
do
        fbname=$(basename "$filename" | cut -d. -f1)
        name1=$(echo ${fbname} | cut -d "_" -f1-6)
        name2=$(echo ${fbname} | cut -d "_" -f1-5)
        mkdir -p ../../Mapping_STAR/${fbname}
        start STAR --runThreadN 32 \
        --genomeDir ../../GenomeIndicesGRCm38.93_STAR \
        --sjdbGTFfile ../../GenomeIndicesGRCm38.93_STAR/Mus_musculus.GRCm38.93.gtf \
        --readFilesIn ${filename} \
        --readFilesCommand zcat \
        --outFileNamePrefix ../../Mapping_STAR/${fbname}/${fbname} \
        --outSAMattributes All  \
        --outSAMstrandField intronMotif \
        --sjdbOverhang 100 \
        --outSAMtype BAM SortedByCoordinate \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --outSJfilterReads All \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --outSAMmultNmax 1 \
        --outSAMattrRGline ID:${name1} SM:${name2}
done
cd ../../

###################################### Merge lanes ############################
cd ./Mapping_STAR
samtools merge -f ./Merged/${name}Aligned.sortedByCoord.out.merged.bam \
               ./${name}_*_L001_R1_001/${name}_*_L001_R1_001Aligned.sortedByCoord.out.bam \
               ./${name}_*_L002_R1_001/${name}_*_L002_R1_001Aligned.sortedByCoord.out.bam 

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



############ run Multiqc ############
cd ./QC/PreFilteringFastQCzip_merged # change directory to PreFilteringFastQCzip_merged
multiqc ./* # do multiQC on report for fastqc on merged preQC files
cd ../PostFilteringFastQCzip_merged # change directory to PreFilteringFastQCzip_merged
multiqc ./* # do multiQC on report for fastqc on merged preQC files
cd ../../ # change directory to reads
done < ../SampleName.txt


rm -rf ./Reads/temp_merged_prefiltering # delete folder temp_merged_prefiltering as that was only for QC check
rm -rf ./Reads/temp_unmerged # delete folder temp_unmerged as that was only for QC check
rm -rf ./Reads/temp_merged_postfiltering # delete folder temp_merged_postfiltering as that was only for QC check