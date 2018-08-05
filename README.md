# RNA-Seq_Mouse
RNA-Seq Analysis Workflow for Mouse Samples

The RNA-Seq reads should be kept in the folder in Reads in the working directory

The workflow consists of the folowing steps:

00_MainNGS.sh: Main script that needs to be run. This script will load all the reequired modules, craete required directories and call the steps for the process. 

01_PreFilteringQC.sh: This script will run FASTqc on original data.


02_FilteringQC.sh: BBDUK suite will be called to trim the reads accoridng to quality of the reads and remove adapters


03_PostFilteringQC.sh: The quality trimmed reads are again subjected to Quality Check using FASTqc


04_Mapping.sh: The reads post trimming are mapped to the reference genome using STAR. Before ampping this script will downlaod the reference genome and annotations as well as create indices for the mapping.


05_MarkDuplicates.sh: The duplicate reads can lead to artifacts so they are removed in this step.


06_FeatureCounts.sh: Now, we need to calulate the reads mapped to genes which is done using FeatureCounts


Rscript DESEQ2.R: The R sctipt is called which used the read counts from previous step and deduces differentially expressed genes