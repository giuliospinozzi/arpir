#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

featCounts <- args[1]
input <- args[2]
output_dir <- args[3]

setwd(output_dir)

# Import data from featureCounts
countdata1 <- read.table(featCounts, header=TRUE, row.names=1,sep = "\t")
lenth_genes <- countdata1[ ,5,drop=F]
countdata <- countdata1[ ,6:ncol(countdata1)]
input_table <- read.csv(input,sep=",")
input_table$sample_name=gsub("-",".",input_table$sample_name)
for (i in 1:nrow(input_table)) {
  colnames(countdata)[grep(input_table$sample_name[i],colnames(countdata))] <- as.character(input_table$sample_name[i])
}

# Convert to matrix
countdata <- as.matrix(countdata)

# experimental design
type <- c()
for (k in 1:ncol(countdata)) {
  type <- c(type,as.character(input_table$Type[grep(colnames(countdata)[k],input_table$sample_name)]))
}
DataGroups <- factor(type)
DataGroups <- relevel(DataGroups,ref=as.character(input_table$Type[1]))

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet
coldata <- data.frame(row.names=colnames(countdata), DataGroups)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~DataGroups)

# PCA
library(ggfortify)
mcols(dds)$basepairs <- as.matrix(lenth_genes)
fpkm <- fpkm(dds)
pca_des <- rbind(fpkm,type=as.character(DataGroups))
png("deseq2-pca.png", w=1000, h=1000, pointsize=30)
autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca_des), colour="type", main="PCA",size=10)+ 
  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=30),
        legend.title=element_blank(),axis.text = element_text(size=30),
        axis.title=element_text(size=30))+geom_text(aes(label=colnames(fpkm)),size=10)
dev.off()

pdf("deseq2-pca.pdf",10,10)
autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca_des), colour="type", main="PCA",size=10)+ 
  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=20),
        legend.title=element_blank(),axis.text = element_text(size=20),
        axis.title=element_text(size=20))+geom_text(aes(label=colnames(fpkm)),size=7)
dev.off()

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Get differential expression results
res <- results(dds)
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(fpkm), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
resdata=resdata[,c(1,3,6:13)]
## Write results
resdata=na.omit(resdata)
resdata[resdata$pvalue==0|resdata$padj==0,c("pvalue","padj")]=0.1e-320
colnames(resdata)[c(2,4)]=c("log2FoldChange","padj")
write.csv(resdata, file="deseq2-diffexpr-results.csv")

## Volcano plot
library(ggrepel)
resdata$color="F"
for (i in 1:nrow(resdata)) {
  if (resdata$padj[i]<0.05) {resdata$color[i]="FDR<0.05"}
  if (abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="|LogFC|>1.5"}
  if (resdata$padj[i]<0.05 & abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="both"}
}
resdata=resdata[order(resdata$padj,decreasing = T),]
png("deseq2-volcanoplot.png", 1200, 1000, pointsize=12)
print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
        geom_point(aes(color = color)) + 
        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
        geom_text_repel(
          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
        theme(plot.title = element_text(hjust = 0.5,face="bold",size=20),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=14),axis.title=element_text(size=14),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.title = element_blank(),legend.text=element_text(size=14),
              legend.background = element_rect(fill="white",
                                               size=0.25, linetype="solid", 
                                               colour ="black")))
dev.off()

pdf("deseq2-volcanoplot.pdf", 15,10)
print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
        geom_point(aes(color = color)) + 
        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
        geom_text_repel(
          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
        theme(plot.title = element_text(hjust = 0.5,face="bold",size=20),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=14),axis.title=element_text(size=14),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.title = element_blank(),legend.text=element_text(size=14),
              legend.background = element_rect(fill="white",
                                               size=0.25, linetype="solid", 
                                               colour ="black")))
dev.off()


# heatmap topgenes
library(genefilter)
library(RColorBrewer)
library(gplots)
topVarGenes <- head( order( rowVars( fpkm ), decreasing=TRUE ), 100 )
png("deseq2-heatmap-topVarGenes.png", w=8, h=9, pointsize=20, res=300, units = "in")
par(cex.main=0.8)
heatmap.2( fpkm[ topVarGenes, ], cexCol=0.5, cexRow=0.3, offsetRow=-0.4, offsetCol=-0.4, 
           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), srtCol=30,
           key.par=list(cex=0.6))
dev.off()

pdf("deseq2-heatmap-topVarGenes.pdf", w=8, h=8)
heatmap.2( fpkm[ topVarGenes, ], cexCol=0.9, cexRow=0.5, offsetRow=-0.4, offsetCol=-0.4, 
           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), srtCol=30)
dev.off()

# Sample distance heatmap
mycols <- brewer.pal(8, "Dark2")[1:length(unique(DataGroups))]
sampleDists <- as.matrix(dist(t(assay(rld))))
png("deseq2-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

pdf("deseq2-heatmap-samples.pdf", w=8, h=8)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
