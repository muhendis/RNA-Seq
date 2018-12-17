source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager",repos = "http://cran.ma.imperial.ac.uk/", dependencies = TRUE,quiet = TRUE)
# Function to Install misisng packages and load required packages
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE,update = TRUE, ask = FALSE,version = "3.8")
    sapply(pkg, require, character.only = TRUE)
}
packages<-c("ggplot2","VennDiagram","EnsDb.Hsapiens.v86","RColorBrewer","pheatmap","DESeq2","varhandle","ggfortify")
check.packages(packages)



setwd("/Users/nsharma/Documents/CRUK-MI_Projects/Alessio/Run_02/new_19_Nov_2018/Mapping_STAR/FeatureCounts_Merged/LnCAP/SetA_Scrbl_MiRNA")
files<-list.files(".",pattern="*Counts.txt$")

#### Read sample data for conditions 
Samples<-read.csv("SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
colnames(Samples)<-c("Sample_ID","condition")
rownames(Samples)<-Samples$Sample_ID

for(i in 1:length(list.files(".",pattern="*Counts.txt$")))
{
  if (!exists("Counts")){
    Counts <- read.table(files[i],header = FALSE, skip=1,stringsAsFactors=FALSE)
    Counts<-Counts[2:nrow(Counts),c(1, ncol(Counts))]
    file.name<-as.character(strsplit(files[i],"Aligned.sortedByCoord.out.merged.DuplicatesRemoved.bam.txt"))
    colnames(Counts)<-c("Gene_ID",file.name)
    next
  }
  if (exists("Counts")){
    temp.Counts <- read.table(files[i],header = FALSE, skip=1,stringsAsFactors=FALSE)
    temp.Counts<-temp.Counts[,c(1, ncol(temp.Counts))]
    file.name<-as.character(strsplit(files[i],"Aligned.sortedByCoord.out.merged.DuplicatesRemoved.bam.txt"))
    colnames(temp.Counts)<-c("Gene_ID",file.name)
    temp.Counts<-temp.Counts[2:nrow(temp.Counts),]
    Counts<-merge(Counts, temp.Counts,by="Gene_ID")
    rm(temp.Counts)
  }
}
row.names(Counts)<-Counts[,1]
Counts<-Counts[,-1]
colnames(Counts)<-c("AC11","AC12","AC13","AC14","AC15","AC16")
Counts.org<-Counts
Counts<-Counts[,c(1,3,5,2,4,6)]
Counts[]<-sapply(Counts[], as.numeric)
Counts<-Counts[which(rowSums(Counts)>0),] # remove the genes with no reads mapped

####### Plot DensityPlot and BoxPlot ########
library(minfi)
library(quantro)
Counts.mat<-as.matrix(Counts)
Counts.mat<-Counts.mat[(rowMeans(Counts.mat)>5),]
Counts.mat<-log2(Counts.mat+1)
matdensity(Counts.mat, groupFactor = Samples$condition, xlab = " ", ylab = "density",main = "Beta Values-LnCAP", brewer.n = 8, brewer.name = "Dark2")
legend('topright', c("Scrbl", "miRNA"), col = c(1, 2), lty = 1, lwd = 3)
matboxplot(Counts.mat, groupFactor = Samples$condition, main = "Beta Values-LnCAP")

##### Deseq2 Object and Normalised Counts ######
dds <- DESeqDataSetFromMatrix(countData = Counts,colData = Samples,design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds)
dds.norm<-estimateSizeFactors(dds)
Counts.normal<-counts(dds.norm, normalized=TRUE)
write.csv(Counts.normal, file = "Counts.normal_DESEQ2.csv")

####### Plot PCA ########
Counts.normal.log<-log2(Counts.normal+1)
Counts.normal.log.t<-as.data.frame(t(Counts.normal.log))
Counts.normal.log.t["type"] = c("scbrl","scbrl","scbrl","miR378a","miR378a","miR378a")
autoplot(prcomp(Counts.normal.log.t[,c(1:(ncol(Counts.normal.log.t)-1))]), 
         data = Counts.normal.log.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4,xlim = c(1, -1), ylim = c(1, -1)) + ggtitle("PCA-LnCAP_31898-Gene")+ theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5))

####### Plot HeatMap ########
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
#TPM.log.isfinite_filtered <- TPM.log.isfinite[apply(TPM.log.isfinite, MARGIN = 1, FUN = function(x) sd(x) != 0),]
Samples2<-as.data.frame(Samples[,2])
names(Samples2)<-"condition"
rownames(Samples2)<-rownames(Samples)
pheatmap(Counts.normal.log, col=colors,show_rownames = FALSE, annotation_col=Samples2, scale="row", main = "LnCAP_Scrbl_MiRNA")

##### Differential analysis using Deseq2  ######
res.miR378aVsscbrl <- results(dds, contrast=c("condition","miR378a","scbrl"), alpha = 0.01) 
write.csv(res.miR378aVsscbrl, file = "res.miR378aVsscbrl.csv")

# Subset significant genes according to Log Fold Change cutoff
res.miR378aVsscbrl.Sig1 <- subset(res.miR378aVsscbrl, res.miR378aVsscbrl$padj <= 0.01 )
res.miR378aVsscbrl.Sig1 <- res.miR378aVsscbrl.Sig1[order(res.miR378aVsscbrl.Sig1$padj),]

res.miR378aVsscbrl.Sig2 <- subset(res.miR378aVsscbrl, res.miR378aVsscbrl$padj <= 0.01 & abs(res.miR378aVsscbrl$log2FoldChange) >=0.75)
res.miR378aVsscbrl.Sig2 <- res.miR378aVsscbrl.Sig2[order(res.miR378aVsscbrl.Sig2$padj),]

res.miR378aVsscbrl.Sig3 <- subset(res.miR378aVsscbrl, res.miR378aVsscbrl$padj <= 0.01 & abs(res.miR378aVsscbrl$log2FoldChange) >=1.5)
res.miR378aVsscbrl.Sig3 <- res.miR378aVsscbrl.Sig3[order(res.miR378aVsscbrl.Sig3$padj),]

# Add annotations (gene names) to the selected genes
convertID<-function(db,ids,key.type,toKey){
  suppressWarnings(x<-mapIds(db, keys=ids, keytype=key.type, column=toKey))
  return(x)
}

res.miR378aVsscbrl.Sig1$Gene_Symbol<-convertID(EnsDb.Hsapiens.v86,row.names(res.miR378aVsscbrl.Sig1), "GENEID","SYMBOL")
write.csv(res.miR378aVsscbrl.Sig1, file = "res.miR378aVsscbrl.padj.01.Sig.csv")

res.miR378aVsscbrl.Sig2$Gene_Symbol<-convertID(EnsDb.Hsapiens.v86,row.names(res.miR378aVsscbrl.Sig2), "GENEID","SYMBOL")
write.csv(res.miR378aVsscbrl.Sig2, file = "res.miR378aVsscbrl.padj.01_lfc.75.Sig.csv")

res.miR378aVsscbrl.Sig3$Gene_Symbol<-convertID(EnsDb.Hsapiens.v86,row.names(res.miR378aVsscbrl.Sig3), "GENEID","SYMBOL")
write.csv(res.miR378aVsscbrl.Sig3, file = "res.miR378aVsscbrl.padj.01_lfc-1.5.Sig.csv")

####### Plot PCA for DEGs ########
Counts.normal.log.subset<-subset(Counts.normal.log,rownames(Counts.normal.log) %in%  rownames(res.miR378aVsscbrl.Sig1))
Counts.normal.log.subset.t<-as.data.frame(t(Counts.normal.log.subset))
Counts.normal.log.subset.t["type"] = c("scbrl","scbrl","scbrl","miR378a","miR378a","miR378a")
autoplot(prcomp(Counts.normal.log.subset.t[,c(1:(ncol(Counts.normal.log.subset.t)-1))]), 
         data = Counts.normal.log.subset.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4,xlim = c(1, -1), ylim = c(1, -1)) + ggtitle("PCA-LnCAP_570-DEGs")+ theme(plot.title = element_text(color="black", size=14, face="bold", hjust=0.5))
