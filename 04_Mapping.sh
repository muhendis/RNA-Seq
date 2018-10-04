#### Create Genome Indice for STAR #### 
cd ./GenomeIndicesGRCm38.93_STAR
wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gzip -d *.gz

STAR --runThreadN 32  \
--readFilesCommand zcat \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
--sjdbGTFfile Mus_musculus.GRCm38.93.gtf

############ Mapping By STAR ############
cd ../QC/QCpassed/
find . -type f -iname '*.fastq.gz' |
while read filename
do
        fbname=$(basename "$filename" | cut -d. -f1)
        name1=$(echo ${fbname} | cut -d "_" -f1-6)
        name2=$(echo ${fbname} | cut -d "_" -f1-5)
        mkdir -p ../../Mapping_STAR/${fbname}
        STAR --runThreadN 32 \
        --genomeDir ../../GenomeIndicesGRCm38.93_STAR \
        --sjdbGTFfile ../../GenomeIndicesGRCm38.93_STAR/Mus_musculus.GRCm38.93.gtf \
        --readFilesIn ${filename} \
        --readFilesCommand zcat \
        --outFileNamePrefix ../../Mapping_STAR/${fbname}/${fbname} \
        --outSAMattributes All  \
        --outSAMstrandField intronMotif \
        --sjdbOverhang 100 \ {related to read length} \
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