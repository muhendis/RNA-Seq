library(varhandle)
library(hash)
library(ggfortify)
library("RColorBrewer")
library("pheatmap")

setwd("/Users/nsharma/Documents/CRUK-MI_Projects/Alessio/Mapping_STAR/FeatureCounts_Duplicates_Removed")

files<-list.files(".",pattern="*.txt$")
genes<-read.table(files[1],header = FALSE, skip=1)
genes<-genes[,c(1, 6)]
colnames(genes)<-c("Gene_ID","Length")
genes<-genes[2:nrow(genes),]
hash()
h<-hash()
k<-as.character(genes[,1])
v<-as.character(genes[,2])
.set(h,keys=c(k),values=c(v))

for(i in 1:length(list.files(".",pattern="*.txt$")))
{
  if (!exists("counts")){
    counts <- read.table(files[i],header = FALSE, skip=1)
    counts<-counts[,c(1, ncol(counts))]
    file.name<-as.character(strsplit(files[i],"Aligned.sortedByCoord.out.merged.bam.txt"))
    colnames(counts)<-c("Gene_ID",file.name)
    counts<-counts[2:nrow(counts),]
    next
  }
  if (exists("counts")){
    temp.counts <- read.table(files[i],header = FALSE, skip=1)
    temp.counts<-temp.counts[,c(1, ncol(temp.counts))]
    file.name<-as.character(strsplit(files[i],"Aligned.sortedByCoord.out.merged.bam.txt"))
    colnames(temp.counts)<-c("Gene_ID",file.name)
    temp.counts<-temp.counts[2:nrow(temp.counts),]
    counts<-merge(counts, temp.counts,by="Gene_ID")
    rm(temp.counts)
  }
}
row.names(counts)<-counts[,1]
counts<-counts[,-1]
colnames(counts)<-c("AC1","AC10","AC2","AC3","AC4","AC5","AC6","AC9")
counts.tmp<-counts
counts<-counts[which(rowSums(unfactor(counts))>0),] # remove the genes with no reads mapped
write.csv(counts,file = "GeneCounts.csv")

###### TPM #######

TPM<-as.data.frame.matrix(counts)

# STEP 1: Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK)
for(i in 1:nrow(TPM))
{
  (gene.name<-rownames(TPM)[i]) # Get name of the gene
  (gene.len<-as.numeric(noquote(h[[ gene.name ]]))/1000)  # Fetch gene length in number of bases converted into kilo bases
  for(j in 1:ncol(TPM))
  {
    TPM[i,j]<-(as.numeric(TPM[i,j]))/gene.len  # divide the number of reads by gene length 
  }
}

# STEP 2: Sum all the reads and divide that by 1000000 to calculate scaling factor
read.sum<-vector()
for(i in 1:ncol(TPM))
{ 
  read.sum[i]<-0
  for(j in 1:nrow(TPM))
  {
    read.sum[i]<-read.sum[i]+(as.numeric(TPM[j,i])) # sum of each column
  }
}
scaling.factor<-read.sum/1000000

# STEP 3: Divide output of step 1 by scaling.factor from step 2
for(i in 1:ncol(TPM))
{ 
  for(j in 1:nrow(TPM))
  {
    TPM[j,i]<-(as.numeric(TPM[j,i]))/scaling.factor[i] # sum of each column
  }
}

TPM<-(TPM[,c(1,3:8,2)])
TPM<-TPM[,c(1,3,5,2,4,6,7,8)]
temp<-TPM
TPM <- as.matrix(as.data.frame(lapply(TPM, as.numeric)))
rownames(TPM)<-rownames(temp)
save(TPM, file = "TPM.RData")

####### TPM calculation complete ########

TPM.log<-log(TPM)
TPM.log.isfinite <- TPM.log[is.finite(rowSums(TPM.log)),]
TPM.log.isfinite.t<-as.data.frame(t(TPM.log.isfinite))
TPM.log.isfinite.t["type"] = c("PC3scbrl", "PC3scbrl", "PC3scbrl", "miR378a", "miR378a", "miR378a", "untreated","untreated")

####### Plot PCA ########

autoplot(prcomp(TPM.log.isfinite.t[,c(1:(ncol(TPM.log.isfinite.t)-1))]), 
         data = TPM.log.isfinite.t, 
         colour = 'type',
         label=TRUE,
         label.size = 4)
dev.off()

####### Plot HeatMap ########

Samples<-read.csv("SampleName.csv",stringsAsFactors=FALSE, header = FALSE)
Samples<-as.data.frame(Samples[,2],row.names = Samples[,1],stringsAsFactors=FALSE)
names(Samples)[1]<-"condition"
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)
pheatmap(TPM.log.isfinite, col=colors,show_rownames = FALSE, annotation_col=Samples, scale="row")
dev.off()

########################################################



