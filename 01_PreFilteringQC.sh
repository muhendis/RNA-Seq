############ Temporary Merging lanes of the samples for QC check  ############
# reads sample names from SampleName.txt and merge the files with output in folder "temp_merged_prefiltering"
cd ./Reads #move to folder with reads which are unmerged and zipped
while read -r line 
do
name="$line"
cat ${name}*.fastq.gz > ./temp_merged_prefiltering/${name}.combined.fastq.gz
done < ../SampleName.txt

############ merged pre filtering  QC check #############
cd ./temp_merged_prefiltering # change directory to temp_merged_prefiltering
find -name '*.fastq.gz' | xargs fastqc -t 32 # pass all merged reads as argument for fastqc
mv *fastqc.zip ../../QC/PreFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv *fastqc.html ../../QC/PreFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
cd ../../QC/PreFilteringFastQCzip_merged # change directory to PreFilteringFastQCzip_merged
multiqc ./* # do multiQC on report for fastqc on merged preQC files
cd ../../ # change directory to reads
rm -rf ./Reads/temp_merged_prefiltering # delete folder temp_merged_prefiltering as that was only for QC check
