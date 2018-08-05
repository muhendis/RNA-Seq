############ Unzip in temp folder ############
cd ./Reads
# unzip all files in folder "Reads" to the the folder "temp_unmerged" while keeping original files intact in reads
find . -type f -iname '*.fastq.gz' | while read filename
do
	fbname=$(basename "$filename")
	s1=`basename "$filename" | cut -d . -f -2`
	gunzip -c -N ${filename}  > ./temp_unmerged/$s1
done
############ run BBDUk for trimming/Filtering on unmerged files ############
# forcetrimleft:remove 11 bases from left end prior to QC, 
# ktrim=r: once a reference kmer is matched in a read, that kmer and all the bases to the right will be trimmed, normal mode for adapter trimming.
# k=13: will match kmers of length 13 or more; ref: to find reference kmers from file adapters.fa
# mink=5:mink allows use of shorter kmers at the ends of the read, so at end of the reads value of k will treated as 5 
# qtrim=rl: Trim read ends to remove bases with quality below trimq from both right and left end
# trimq=10: Regions with average quality BELOW this will be trimmed; hdist=1: hdist" means "hamming distance"; this allows one mismatch
# minlength=20: Reads shorter than this (20 bases in our case) after trimming will be discarded
# stats: Write statistics about which contamininants were detected overwrite=t: Grant permission to overwrite files.
cd ./temp_unmerged  # change directory to "temp_unmerged" to run QC check using BBDUK 
find . -type f -iname '*.fastq' | while read filename
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
done # the QC passed and failed reads are saved in folder /QC/QCpassed /QC/QCfailed respectively
cd ../../ # change directory to reads
rm -rf ./Reads/temp_unmerged # delete folder temp_unmerged as that was only for QC check
