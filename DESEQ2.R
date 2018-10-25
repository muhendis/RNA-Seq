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
packages<-c("gplots","VennDiagram","EnsDb.Mmusculus.v79","RColorBrewer","pheatmap","DESeq2","varhandle","ggfortify")
check.packages(packages)

setwd("/Users/nsharma/Documents/CRUK-MI_Projects/Alessio/Run_01/Mapping_STAR/FeatureCounts_Deduped_Paired")
files<-list.files(".",pattern="*.txt$")

#### Read sample data for conditions 
Samples<-read.csv("SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
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
colnames(Counts)<-c("AC1","AC10","AC2","AC3","AC4","AC5","AC6","AC9")
Counts[]<-sapply(Counts[], as.numeric)
Counts<-Counts[,c(c(1,3:8,2))]
Counts<-Counts[,c(1,3,5,2,4,6,7,8)]
Counts.org<-Counts
Counts<-Counts[which(rowSums(Counts)>0),] # remove the genes with no reads mapped
Samples.org<-Samples
Samples$ID<-rownames(Samples)

##################### Differential expression analysis w.r.t "untreated" ##################### 
dds <- DESeqDataSetFromMatrix(countData = Counts,colData = Samples,design = ~ condition)
dds <- DESeq(dds)

### calculate normalised Counts ###
dds.norm<-estimateSizeFactors(dds)
Counts.normal<-counts(dds.norm, normalized=TRUE)
Counts.normal.org<-Counts.normal
write.csv(Counts.normal.org,"Normalised_Gene_Counts_Deseq2.csv")
##################### PLotting using Normalised data from Deseq2  for all samples ############################
Counts.normal.log<-log(Counts.normal)
Counts.normal.log.isfinite <- Counts.normal.log[is.finite(rowSums(Counts.normal.log)),]
Counts.normal.log.isfinite.t<-as.data.frame(t(Counts.normal.log.isfinite))
Counts.normal.log.isfinite.t["type"] = c("scbrl", "scbrl", "scbrl", "miR378a", "miR378a", "miR378a", "untreated","untreated")
####### Plot PCA ########
autoplot(prcomp(Counts.normal.log.isfinite.t[,c(1:(ncol(Counts.normal.log.isfinite.t)-1))]), 
         data = Counts.normal.log.isfinite.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4)
####### Plot HeatMap ########
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(Counts.normal.log.isfinite, col=colors,show_rownames = FALSE, annotation_col=Samples, scale="row")

#######################################################

res.miR378a.untreated <- results(dds, contrast=c("condition","miR378a","untreated"), alpha = 0.05)  ##### miR378a vs untreated ####
res.scbrl.untreated <- results(dds, contrast=c("condition","scbrl","untreated"), alpha = 0.05)   ##### scbrl vs untreated ####
res.miR378a.scbrl <- results(dds, contrast=c("condition","miR378a","scbrl"), alpha = 0.05)  ##### miR378a vs scbrl ####

# Create columns for inforamation on significant genes later used in volcano plot
res.miR378a.untreated<-as.data.frame(dplyr::mutate(as.data.frame(res.miR378a.untreated), sig=ifelse(res.miR378a.untreated$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.miR378a.untreated))
res.scbrl.untreated<-as.data.frame(dplyr::mutate(as.data.frame(res.scbrl.untreated), sig=ifelse(res.scbrl.untreated$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.scbrl.untreated))
res.miR378a.scbrl<-as.data.frame(dplyr::mutate(as.data.frame(res.miR378a.scbrl), sig=ifelse(res.miR378a.scbrl$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res.miR378a.scbrl))
write.csv(res.miR378a.untreated, file = "res.miR378a.untreated.csv")
write.csv(res.scbrl.untreated, file = "res.scbrl.untreated.csv")
write.csv(res.miR378a.scbrl, file = "res.miR378a.scbrl.csv")

# Subset significant genes according to Log Fold Change cutoff
res.miR378a.untreated.Sig <- subset(res.miR378a.untreated, res.miR378a.untreated$padj < 0.05 & abs(res.miR378a.untreated$log2FoldChange) >=1)
res.miR378a.untreated.Sig <- res.miR378a.untreated.Sig[order(res.miR378a.untreated.Sig$padj),]
res.scbrl.untreated.Sig <- subset(res.scbrl.untreated, res.scbrl.untreated$padj < 0.05 & abs(res.scbrl.untreated$log2FoldChange) >=1)
res.scbrl.untreated.Sig <- res.scbrl.untreated.Sig[order(res.scbrl.untreated.Sig$padj),]
res.miR378a.scbrl.Sig <- subset(res.miR378a.scbrl, res.miR378a.scbrl$padj < 0.05 & abs(res.miR378a.scbrl$log2FoldChange) >=1)
res.miR378a.scbrl.Sig <- res.miR378a.scbrl.Sig[order(res.miR378a.scbrl.Sig$padj),]
write.csv(res.miR378a.untreated.Sig, file = "res.miR378a.untreated.Sig.csv")
write.csv(res.scbrl.untreated.Sig, file = "res.scbrl.untreated.Sig.csv")
write.csv(res.miR378a.scbrl.Sig, file = "res.miR378a.scbrl.Sig.csv")

##################### PLotting using Normalised data from Deseq2 without untreated ############################
Counts.normal<-Counts.normal[,c(1:6)]
Samples<-Samples[1:6,]
Counts.normal.log<-log(Counts.normal)
Counts.normal.log.isfinite <- Counts.normal.log[is.finite(rowSums(Counts.normal.log)),]
Counts.normal.log.isfinite.t<-as.data.frame(t(Counts.normal.log.isfinite))
Counts.normal.log.isfinite.t["type"] = c("scbrl", "scbrl", "scbrl", "miR378a", "miR378a", "miR378a")
####### Plot PCA ########
autoplot(prcomp(Counts.normal.log.isfinite.t[,c(1:(ncol(Counts.normal.log.isfinite.t)-1))]), 
         data = Counts.normal.log.isfinite.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4)
####### Plot HeatMap ########
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(Counts.normal.log.isfinite, col=colors,show_rownames = FALSE, annotation_col=Samples, scale="row")

############### for the list of differentialy expressed genes between miR378a Vs scbrl plot PCA and HeatMap 
Counts.normal<-Counts.normal.org
Samples<-Samples.org
sig.genes.matrix<- res.miR378a.scbrl.Sig[,1:6]
Counts.normal<- Counts.normal[which(rownames(Counts.normal) %in% rownames(sig.genes.matrix)),]
Counts.normal.log<-log(Counts.normal)
Counts.normal.log.isfinite <- Counts.normal.log[is.finite(rowSums(Counts.normal.log)),]
Counts.normal.log.isfinite.t<-as.data.frame(t(Counts.normal.log.isfinite))
Counts.normal.log.isfinite.t["type"] = c("scbrl", "scbrl", "scbrl", "miR378a", "miR378a", "miR378a", "untreated","untreated")
##################### Differential expression analysis w.r.t "untreated" ##################### 
####### Plot PCA ########
autoplot(prcomp(Counts.normal.log.isfinite.t[,c(1:(ncol(Counts.normal.log.isfinite.t)-1))]), 
         data = Counts.normal.log.isfinite.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4)
####### Plot HeatMap ########
Samples<-read.csv("SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
Samples<-as.data.frame(Samples[,2],row.names = Samples[,1],stringsAsFactors=FALSE)
names(Samples)[1]<-"condition"
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(Counts.normal.log.isfinite, col=colors,show_rownames = T, fontsize_row=5, annotation_col=Samples, scale="row")

##################### PLotting using Normalised data from Deseq2 without untreated ############################
Counts.normal<-Counts.normal[,c(1:6)]
Counts.normal.log<-log(Counts.normal)
Counts.normal.log.isfinite <- Counts.normal.log[is.finite(rowSums(Counts.normal.log)),]
Counts.normal.log.isfinite.t<-as.data.frame(t(Counts.normal.log.isfinite))
Counts.normal.log.isfinite.t["type"] = c("scbrl", "scbrl", "scbrl", "miR378a", "miR378a", "miR378a")
####### Plot PCA ########
autoplot(prcomp(Counts.normal.log.isfinite.t[,c(1:(ncol(Counts.normal.log.isfinite.t)-1))]), 
         data = Counts.normal.log.isfinite.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4)
####### Plot HeatMap ########
Samples<-read.csv("SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
Samples<-as.data.frame(Samples[,2],row.names = Samples[,1],stringsAsFactors=FALSE)
names(Samples)[1]<-"condition"
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(Counts.normal.log.isfinite, col=colors,show_rownames = T, fontsize_row=5,annotation_col=Samples, scale="row",width = 12, height = 6)






######################################## Moule to annotate mouse genome ##############################
# Add annotations (gene names) to the selected genes
#convertID<-function(db,ids,key.type,toKey){
#  suppressWarnings(x<-mapIds(db, keys=ids, keytype=key.type, column=toKey))
#  return(x)
#}
#res.RT.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.RT.Sig), "GENEID","SYMBOL")
#res.Combo.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.Combo.Sig), "GENEID","SYMBOL")
#res.CD40.Sig$Gene_ID<-convertID(EnsDb.Mmusculus.v79,row.names(res.CD40.Sig), "GENEID","SYMBOL")
####################################################################################################


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

