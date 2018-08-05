############ merged  QC passed reaads #############
cd ./QC/QCpassed # change directory to QCpassed
# reads sample names from SampleName.txt and merge the files with output in folder "temp_merged_prefiltering"
while read -r line 
do
name="$line"
echo $name; echo "..."
cat ${name}*.fastq.gz > ../../Reads/temp_merged_postfiltering/${name}.QCpassed.combined.fastq.gz
done < ../../SampleName.txt

########### FastQC for merged QC passed reads ############
cd ../../Reads/temp_merged_postfiltering # change directory to temp_merged_prefiltering
find -name '*.fastq.gz' | xargs fastqc -t 32 # pass all merged reads as argument for fastqc
mv *fastqc.zip ../../QC/PostFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv *fastqc.html ../../QC/PostFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
cd ../../QC/PostFilteringFastQCzip_merged # change directory to PreFilteringFastQCzip_merged
multiqc ./* # do multiQC on report for fastqc on merged preQC files
cd ../../ # change directory to reads
rm -rf ./Reads/temp_merged_postfiltering # delete folder temp_merged_postfiltering as that was only for QC check