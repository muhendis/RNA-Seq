#!/bin/bash
# Filter
#PBS -N Filter
#PBS -l nodes=1:large,walltime=5:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogFilter
#PBS -e /scratch/wsspaces/nsharma-test-0/LogFilter
#PBS -j oe
#PBS -V

############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############

module load apps/bbmap/36.20
module load apps/fastqc/0.11.3/linux-x86_64
module load apps/multiqc/1.4

############ Create Directories #############

mkdir -p ./QC
mkdir -p ./QC/PostFilteringFastQCzip
mkdir -p ./QC/PostFilteringFastQChtml
mkdir -p ./QC/QCpassed
mkdir -p ./QC/QCfailed
mkdir -p ./Reads/temp_unmerged

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
############ run BBDUk for trimming/Filtering on unmerged files ############
cd ./temp_unmerged  # change directory to "temp_unmerged" to run QC check using BBDUK
 for i in {4..5}
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
   ktrim=r \
   k=13 \
   mink=5 \
   qtrim=rl \
   trimq=10 \
   hdist=1 \
   minlength=20 \
   stats=../../QC/bbduk_${name}_L00${i}.stats.txt
done
done < ../SampleName.txt
cd ../
rm -rf temp_unmerged
########### FastQC for merged QC passed reads ############
cd ../QC/QCpassed
find -name '*.fastq.gz' | xargs fastqc -t 32  # pass all reads as argument for fastqc
mv *fastqc.zip ../PostFilteringFastQCzip/ 
mv *fastqc.html ../PostFilteringFastQChtml/ 
cd ../PostFilteringFastQCzip
multiqc ./* 

