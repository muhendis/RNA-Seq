# Bulk RNA-Seq Analysis Workflow

The RNA-Seq reads should be kept in the folder in Reads in the working directory

The workflow consists of the folowing steps:

Step 1: 01_PreFilteringQC.sh: To check the quality of data

Step 2:  a) If the quality of data is not good accorind to Step 1 then 
02_Filtering.sh.sh: BBDUK suite will be called to trim the reads accoridng to quality of the reads and remove adapters and the quality trimmed reads is analysed using FASTqc

b) If the quality of data is  good accorind to Step 1 then
03_StarCreateGenomeIndices.sh: to create indices for reference genome and 
04_Mapping.sh: Map reads to the reference genome using STAR

Step 3:
05_FeatureCounts_Original_Paired.sh: Calculate feature using feature counts.

Step 4:
06_DESEQ2.R: The R script is called which used the read counts from previous step and deduces differentially expressed genes
