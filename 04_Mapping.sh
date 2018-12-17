#!/bin/bash
# Mapping
#PBS -N Mapping
#PBS -l nodes=1:large,walltime=12:00:00
#PBS -m bea
#PBS -M nitin.sharma@cruk.manchester.ac.uk
#PBS -o /scratch/wsspaces/nsharma-test-0/LogMap
#PBS -e /scratch/wsspaces/nsharma-test-0/LogMap
#PBS -j oe
#PBS -V

############ Change Directory ############

cd /scratch/wsspaces/nsharma-test-0/

############ Load modules ############

module load apps/star/2.5.1b/gcc-5.1.0
module load apps/samtools/1.3.1/gcc-4.4.7
module load apps/multiqc/1.4

############ Create Directories #############

mkdir -p ./Mapping_STAR
mkdir -p ./Mapping_STAR/multiqc
mkdir -p ./tmp

############## Main Script #######
cd ./Reads
while read -r line
do
name="$line"
############ Mapping By STAR ############
for (( i=4; i<=5; i++ )){
STAR --genomeDir ../STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94 \
     --sjdbGTFfile ../STAR_GenomeIndices_sjdbOverhang_100_Ensembl_GRCh38.94/Homo_sapiens.GRCh38.94.gtf \
     --readFilesIn ${name}_L00${i}_R1_001.fastq.gz ${name}_L00${i}_R2_001.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ../Mapping_STAR/${name}_L00${i} \
     --sjdbOverhang 100 \
     --outSAMtype BAM SortedByCoordinate \
     --twopassMode Basic \
     --quantMode GeneCounts \
     --outFilterMultimapNmax 1 \
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
 }
done < ../SampleName.txt

cd ../Mapping_STAR
find . -type f -iname "[a-zA-Z0-9]*Log.final.out" | xargs cp -t ./multiqc
find . -type f -iname "[a-zA-Z0-9]*ReadsPerGene.out.tab" | xargs cp -t ./multiqc
cd ./multiqc
multiqc ./* 