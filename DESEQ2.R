source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager",repos = "http://cran.ma.imperial.ac.uk/", dependencies = TRUE,quiet = TRUE)
BiocManager::install(update = TRUE, ask = FALSE) #    Installed packages can be updated to their current version
# Function to Install misisng packages and load required packages
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        BiocManager::install(new.pkg, dependencies = TRUE,update = TRUE, ask = FALSE)
    sapply(pkg, require, character.only = TRUE)
}
packages<-c("gplots","VennDiagram","EnsDb.Mmusculus.v79","RColorBrewer","pheatmap")
check.packages(packages)

setwd("./FeatureCounts")
files<-list.files(".",pattern="*.txt$")

#### Read sample data for conditions 
Samples<-read.csv("../SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
Samples<-as.data.frame(Samples[,2],row.names = Samples[,1],stringsAsFactors=FALSE)
names(Samples)[1]<-"condition"

for(i in 1:length(list.files(".",pattern="*.txt$")))
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
colnames(Counts)<-c("TI2_01","TI2_02","TI2_03","TI2_04","TI2_05","TI2_06","TI2_07","TI2_08","TI2_09","TI2_10","TI2_11","TI2_12","TI2_13","TI2_14","TI2_15","TI2_16")
Counts[]<-sapply(Counts[], as.numeric)
Counts<-Counts[,c(1:3,5:10,12)]
Counts<-Counts[which(rowSums(Counts)>0),] # remove the genes with no reads mapped

Samples.org<-Samples
Counts.org<-Counts

Samples$ID<-rownames(Samples)
Samples<-as.data.frame(Samples[-c(1,7),])
Counts<-Counts[,-c(1,7)]

library("DESeq2")
##### Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = Counts,colData = Samples,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Untreated")
dds <- DESeq(dds)
#resultsNames(dds)

### calculate normalised Counts ###
dds.norm<-estimateSizeFactors(dds)
Counts.normal<-counts(dds.norm, normalized=TRUE)

res.CD40 <- results(dds, contrast=c("condition","CD40","Untreated"), alpha = 0.05)  ##### CD40 vs Control ####
res.RT <- results(dds, contrast=c("condition","RT","Untreated"), alpha = 0.05)   ##### RT vs Control ####
res.Combo <- results(dds, contrast=c("condition","Combination","Untreated"), alpha = 0.05)##### CD40+RT vs Control ####

# Create columns for inforamation on significant genes later used in volcano plot
res.Combo<-as.data.frame(dplyr::mutate(as.data.frame(res.Combo), sig=ifelse(res.Combo$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.Combo))
res.RT<-as.data.frame(dplyr::mutate(as.data.frame(res.RT), sig=ifelse(res.RT$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.RT))
res.CD40<-as.data.frame(dplyr::mutate(as.data.frame(res.CD40), sig=ifelse(res.CD40$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.CD40))

# Subset significant genes according to Log Fold Change cutoff
res.CD40.Sig <- subset(res.CD40, res.CD40$padj < 0.05 & abs(res.CD40$log2FoldChange) >=1); res.CD40.Sig <- res.CD40.Sig[order(res.CD40.Sig$padj),]
res.Combo.Sig <- subset(res.Combo, res.Combo$padj < 0.05 & abs(res.Combo$log2FoldChange) >=1); res.Combo.Sig <- res.Combo.Sig[order(res.Combo.Sig$padj),]
res.RT.Sig <- subset(res.RT, res.RT$padj < 0.05 & abs(res.RT$log2FoldChange) >=1); res.RT.Sig <- res.RT.Sig[order(res.RT.Sig$padj),]

# Add annotations (gene names) to the selected genes
convertID<-function(db,ids,key.type,toKey){
  suppressWarnings(x<-mapIds(db, keys=ids, keytype=key.type, column=toKey))
  return(x)
}

res.RT.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.RT.Sig), "GENEID","SYMBOL")
res.Combo.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.Combo.Sig), "GENEID","SYMBOL")
res.CD40.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.CD40.Sig), "GENEID","SYMBOL")
write.csv(res.CD40.Sig, file = "res.CD40.Sig.csv")
write.csv(res.Combo.Sig, file = "res.Combo.Sig.csv")
write.csv(res.RT.Sig, file = "res.RT.Sig.csv")


## significant genes combined
sig.genes<-unique(c(row.names(res.CD40.Sig),row.names(res.RT.Sig),row.names(res.Combo.Sig)))
## normalised counts for the significant genes
sig.genes.counts<-Counts.normal[which(rownames(Counts.normal) %in% sig.genes),]
### heatmap using the normalised counts 
x<-log2(sig.genes.counts)
x[is.infinite(x)] <- 0
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(x, col=colors,show_rownames = FALSE, annotation_col=Samples, scale="row")
dev.off()

res.Combo.Sig.down<-subset(res.Combo.Sig,log2FoldChange<0,select=Gene_ID)
res.Combo.Sig.up<-subset(res.Combo.Sig,log2FoldChange>0,select=Gene_ID)
res.RT.Sig.up<-subset(res.RT.Sig,log2FoldChange>0,select=Gene_ID)
res.RT.Sig.down<-subset(res.RT.Sig,log2FoldChange<0,select=Gene_ID)
res.CD40.Sig.down<-subset(res.CD40.Sig,log2FoldChange<0,select=Gene_ID)
res.CD40.Sig.up<-subset(res.CD40.Sig,log2FoldChange>0,select=Gene_ID)
venn.plot <- venn.diagram(list(row.names(res.CD40.Sig.up), row.names(res.RT.Sig.up), row.names(res.Combo.Sig.up)), NULL, main="Up regulated genes", fill=c("red", "blue","yellow"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CD40", "RT","Combo"))
grid.draw(venn.plot)
dev.off()
venn.plot <- venn.diagram(list(row.names(res.CD40.Sig.down), row.names(res.RT.Sig.down), row.names(res.Combo.Sig.down)), NULL, main="Down-Regulated genes", fill=c("red", "blue","yellow"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CD40", "RT","Combo"))
grid.draw(venn.plot)
dev.off()
venn.plot <- venn.diagram(list(row.names(res.CD40.Sig), row.names(res.RT.Sig), row.names(res.Combo.Sig)), NULL, main="All Significant genes", fill=c("red", "blue","yellow"), alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, category.names=c("CD40", "RT","Combo"))
grid.draw(venn.plot)
dev.off()

#### Log fold change shrinkage for visualization and ranking
res.CD40.LFC <- lfcShrink(dds, coef="condition_CD40_vs_Untreated", type="apeglm")
res.RT.LFC <- lfcShrink(dds, coef="condition_RT_vs_Untreated", type="apeglm")
res.Combo.LFC <- lfcShrink(dds, coef="condition_Combo_vs_Untreated", type="apeglm")


### MA-plot
plotMA(res.Combo.LFC, ylim=c(-2,2))
plotMA(res.CD40.LFC, ylim=c(-2,2))
plotMA(res.RT.LFC, ylim=c(-2,2))

#### Check the expression levels of the most differentially expressed gene
most.sign.gene <- rownames(res.Combo.Sig)[1]
#gn.most.sign<-"ENSMUSG00000000001"
most.sign.gene.expression <- counts(dds.norm, normalized=T)[most.sign.gene,]
barplot(most.sign.gene.expression, main=most.sign.gene, las=2, cex.names=0.5, col = c("Red","Red","Red","Blue","Blue","Blue","Blue","Green","Green","Green","Orange","Orange","Orange","Orange"))
#plotCounts(dds, gene=which.min(resCD40$padj), intgroup="condition")
#plotCounts(dds, gene="ENSMUSG00000034438", intgroup="condition")

##### Volcano Plot
res.Combo <- res.Combo[order(res.Combo$padj),] # order results table by the smallest adjusted p value
results.Combo<-as.data.frame(dplyr::mutate(as.data.frame(res.Combo), sig=ifelse(res.Combo$padj<0.01, "FDR<0.01", "Not Sig")), row.names=rownames(res.Combo))
results<-results.Combo
DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(2) & results$padj < 0.01),]
p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")
p + ggrepel::geom_text_repel(data=results[1:5, ], ggplot2::aes(label=rownames(results[1:5, ])))

##################################################################################################

