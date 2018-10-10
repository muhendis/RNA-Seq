#!/bin/bash
# NGS1
#PBS -N NGS1
#PBS -l nodes=1:ppn=16,mem=90154756kb,walltime=36:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-Alessio_ACseries-0/Log1
#PBS -e /scratch/wsspaces/nsharma-Alessio_ACseries-0/Log1
#PBS -j oe
#PBS -V
##PBS -A tartheonc

############ Change Directory ############

cd /scratch/wsspaces/nsharma-Alessio_ACseries-0/

############ Load modules ############

module load apps/bbmap/36.20
module load apps/fastqc/0.11.3/linux-x86_64
module load apps/multiqc/1.4
module load apps/star/2.5.1b/gcc-5.1.0
module load apps/picardtools/1.96/noarch
module load apps/samtools/1.3.1/gcc-4.4.7
module load apps/subread/1.5.0-p3/gcc-4.4.7

############ Create Directories #############

mkdir -p ./QC
mkdir -p ./QC/PreFilteringFastQCzip_merged
mkdir -p ./QC/PreFilteringFastQChtml_merged
mkdir -p ./QC/PostFilteringFastQCzip_merged
mkdir -p ./QC/PostFilteringFastQChtml_merged
mkdir -p ./QC/QCpassed
mkdir -p ./QC/QCfailed
mkdir -p ./Reads/temp_unmerged
mkdir -p ./Reads/temp_merged_prefiltering
mkdir -p ./Reads/temp_merged_postfiltering
mkdir -p ./Mapping_STAR
mkdir -p ./Mapping_STAR/Merged
mkdir -p ./Mapping_STAR/Duplicates_Removed
mkdir -p ./Mapping_STAR/multiqc
mkdir -p ./tmp

############## Main Script #######
cd ./Reads
while read -r line
do
name="$line"
############ Unzip in temp folder ############
# unzip all files in folder "Reads" to the the folder "temp_unmerged" while keeping original files intact in reads
find . -type f -iname "${name}*.fastq.gz" | while read filename
do
    fbname=$(basename "$filename")
    s1=`basename "$filename" | cut -d . -f -2`
    gunzip -c -N ${filename}  > ./temp_unmerged/$s1
done
############ Merge lanes of the samples for QC check  ############
cat ${name}*R1*.fastq.gz > ./temp_merged_prefiltering/${name}.R1.combined.fastq.gz
cat ${name}*R2*.fastq.gz > ./temp_merged_prefiltering/${name}.R2.combined.fastq.gz
############ merged pre filtering  QC check #############
cd ./temp_merged_prefiltering
find -name "${name}*.combined.fastq.gz" | xargs fastqc -t 16
mv ${name}*fastqc.zip ../../QC/PreFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv ${name}*fastqc.html ../../QC/PreFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
############ run BBDUk for trimming/Filtering on unmerged files ############
cd ../temp_unmerged  # change directory to "temp_unmerged" to run QC check using BBDUK
 for i in {1..6}
  do
   bbduk.sh \
   in=${name}_L00${i}_R1_001.fastq \
   in2=${name}_L00${i}_R2_001.fastq \
   out=../../QC/QCpassed/${name}_L00${i}_R1_001.fastq.gz \
   out2=../../QC/QCpassed/${name}_L00${i}_R2_001.fastq.gz \
   outm=../../QC/QCfailed/${name}_L00${i}_R1_001.fastq.gz \
   outm2=../../QC/QCfailed/${name}_L00${i}_R2_001.fastq.gz \
   ref=/apps/modules/pkg/apps/bbmap/36.20/resources/adapters.fa \
   overwrite=t \
   forcetrimleft=11 \
   ktrim=r \
   k=13 \
   mink=5 \
   qtrim=rl \
   trimq=10 \
   hdist=1 \
   minlength=20 \
   stats=../../QC/bbduk_${name}_L00${i}.stats.txt
done
############ merged  QC passed reaads #############
cd ../../QC/QCpassed # change directory to QCpassed
# reads sample names from SampleName.txt and merge the files with output in folder "temp_merged_prefiltering"
cat ${name}*R1*.fastq.gz > ../../Reads/temp_merged_postfiltering/${name}.R1.QCpassed.combined.fastq.gz
cat ${name}*R2*.fastq.gz > ../../Reads/temp_merged_postfiltering/${name}.R2.QCpassed.combined.fastq.gz
########### FastQC for merged QC passed reads ############
cd ../../Reads/temp_merged_postfiltering # change directory to temp_merged_prefiltering
find -name "${name}*.QCpassed.combined.fastq.gz" | xargs fastqc -t 16 # pass all merged reads as argument for fastqc
mv ${name}*fastqc.zip ../../QC/PostFilteringFastQCzip_merged/ # move fastqc zip file to the folder PreFilteringFastQCzip_merged
mv ${name}*fastqc.html ../../QC/PostFilteringFastQChtml_merged/ # move fastqc html file to the folder PreFilteringFastQChtml_merged
############ Mapping By STAR ############
mkdir -p ../../Mapping_STAR/${name}
STAR --runThreadN 16 \
      --genomeDir ../../GenomeIndicesSTAR_hg38 \
     --sjdbGTFfile ../../GenomeIndicesSTAR_hg38/hg38.gtf \
     --readFilesIn ${name}.R1.QCpassed.combined.fastq.gz ${name}.R2.QCpassed.combined.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ../../Mapping_STAR/${name} \
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
     --outFilterType BySJout \
     --outWigType bedGraph \
     --outWigStrand Stranded \
     --outWigNorm RPM \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --alignIntronMin 20 \
     --alignIntronMax 1000000 \
     --alignMatesGapMax 1000000 \
     --outFilterMismatchNoverReadLmax 0.04
cd ../../Mapping_STAR
cp *Log.final.out ./multiqc
cp *ReadsPerGene.out.tab ./multiqc
############################ QC_Picardtools ############################
find . -type f -iname '${name}*Aligned.sortedByCoord.out.bam' |
while read filename
do
          fbname=$(basename "$filename" | cut -d. -f1)
          postfix=".sortedByCoord.out.Deduped.bam"
          java -Xmx64g \
          -Djava.io.tmpdir=../tmp \
          -XX:ParallelGCThreads=8 \
          -jar /home/nsharma/anaconda3/share/picard-2.18.9-0/picard.jar \
          MarkDuplicates \
          REMOVE_DUPLICATES=true \
          INPUT=$filename \
          OUTPUT=./Duplicates_Removed/${fbname}${postfix} \
          METRICS_FILE=./Duplicates_Removed/${fbname}_mark_dups_metrics.txt
done
cd ../
done < ../SampleName_1.txt
