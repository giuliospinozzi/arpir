#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

cuffdiff <- args[1]
input <- args[2]
output_dir <- args[3]

setwd(output_dir)

library(cummeRbund)
cuff_data <- readCufflinks(cuffdiff,rebuild = T)

# PCA
library(ggfortify)
input_table <- read.csv(input,sep=",")
fpkm <- repFpkmMatrix(genes(cuff_data))
#change sample names
sam_names=read.table(paste0(cuffdiff,"/read_groups.info"),header = T)
sam_names$name1=paste0(sam_names$condition,"_",sam_names$replicate_num)
for (i in 1:nrow(sam_names)) {
  sam_names$name2[i]=gsub("^.*?Quantification_and_DEA/","",sam_names$file[i])
  sam_names$name2[i]=gsub("/cuffquant/abundances.cxb$","",sam_names$name2[i])
  colnames(fpkm)[grep(sam_names$name1[i],colnames(fpkm))]=sam_names$name2[i]
}
DataGroups <- c()
for (i in 1:ncol(fpkm)) {
  DataGroups=append(DataGroups,as.character(input_table[input_table$sample_name==colnames(fpkm)[i],"Type"]))
}
pca <- rbind(fpkm,type=as.character(DataGroups))
png("cummeRbund-pca.png", w=1000, h=1000, pointsize=30)
autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca), colour="type", main="PCA",size=10)+ 
  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=30),
        legend.title=element_blank(),axis.text = element_text(size=30),
        axis.title=element_text(size=30))+geom_text(aes(label=colnames(fpkm)),size=10)
dev.off()

pdf("cummeRbund-pca.pdf", 10,10)
autoplot(prcomp(log2((t(fpkm))+1)),data=t(pca), colour="type", main="PCA",size=10)+ 
  theme(plot.title = element_text(face="bold",hjust=0.5,size=50),legend.text=element_text(size=20),
        legend.title=element_blank(),axis.text = element_text(size=20),
        axis.title=element_text(size=20))+geom_text(aes(label=colnames(fpkm)),size=7)
dev.off()

# #Retrive significant gene IDs (XLOC) with a pre-specified alpha
diffGeneIDs <- diffData(genes(cuff_data))
#Use returned identifiers to create a CuffGeneSet object with all relevant info for given genes
diffGenes <- cummeRbund::getGenes(cuff_data,diffGeneIDs$gene_id)

#gene_short_name values (and corresponding XLOC_* values) can be retrieved from the CuffGeneSet by using:
names <- featureNames(diffGenes)
row.names(names) <- names$tracking_id
diffGenesNames <- as.matrix(names)
diffGenesNames <- diffGenesNames[,-1]

# get the data for the significant genes
diffGenesData <- diffData(diffGenes)
row.names(diffGenesData) <- diffGenesData$gene_id
diffGenesData <- diffGenesData[,-1]

# merge the two matrices by row names
diffGenesOutput <- merge(diffGenesNames,diffGenesData,by="row.names")
rownames(diffGenesOutput) <- diffGenesOutput$Row.names
diffGenesOutput <- diffGenesOutput[,-1]
diffGenesOutput1 <- merge(diffGenesOutput,fpkm,by="row.names")
diffGenesOutput1=diffGenesOutput1[,c(1,2,8,10,13:ncol(diffGenesOutput1))]
colnames(diffGenesOutput1)[2:4] <- c("Gene","log2FoldChange","padj")
diffGenesOutput1 <- diffGenesOutput1[order(diffGenesOutput1$padj), ]
diffGenesOutput1=na.omit(diffGenesOutput1)
diffGenesOutput1[diffGenesOutput1$padj==0,c("padj")]=0.1e-320
write.csv(diffGenesOutput1[,2:ncol(diffGenesOutput1)], file="cummeRbund-diffexpr-results.csv")

# volcano plot
library(ggrepel)
resdata=diffGenesOutput1
resdata$color="F"
for (i in 1:nrow(resdata)) {
  if (resdata$padj[i]<0.05) {resdata$color[i]="FDR<0.05"}
  if (abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="|LogFC|>1.5"}
  if (resdata$padj[i]<0.05 & abs(resdata$log2FoldChange[i])>1.5) {resdata$color[i]="both"}
}
resdata=resdata[order(resdata$padj,decreasing = T),]
png("cummeRbund-volcanoplot.png", 1200, 1000, pointsize=12)
print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
        geom_point(aes(color = color)) + 
        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
        geom_text_repel(
          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
        theme(plot.title = element_text(hjust = 0.5,face="bold",size=30),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=14),axis.title=element_text(size=14),
              legend.key = element_rect(colour = NA, fill = NA),
              legend.title = element_blank(),legend.text=element_text(size=14),
              legend.background = element_rect(fill="white",
                                               size=0.25, linetype="solid", 
                                               colour ="black")))
dev.off()

pdf("cummeRbund-volcanoplot.pdf", 15,10)
print(ggplot(resdata, aes(log2FoldChange, -log10(padj)))+
        geom_point(aes(color = color)) + 
        scale_color_manual(values = c("F"="black","FDR<0.05"="red","|LogFC|>1.5"="orange","both"="green"),
                           breaks=c("F","FDR<0.05", "|LogFC|>1.5","both"), name="",
                           labels=c("F","FDR<0.05", "|LogFC|>1.5","both"),
                           limits=c("FDR<0.05", "|LogFC|>1.5","both")) +
        geom_text_repel(
          data = subset(resdata, abs(log2FoldChange)>3 & padj <.05),
          aes(label = Gene), size=2, segment.size=0.2) + ggtitle("Volcano Plot") +
        theme(plot.title = element_text(hjust = 0.5,face="bold",size=30),
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
fpkm <- as.matrix(fpkm)
topVarGenes <- head( order( rowVars( fpkm ), decreasing=TRUE ), 100 )
rownames(diffGenesOutput1)=diffGenesOutput1$Row.names
lab <- as.character(diffGenesOutput1[rownames(fpkm[ topVarGenes, ]),"Gene"])
png("cummeRbund-heatmap-topVarGenes.png", w=8, h=9, pointsize=20, res=300, units = "in")
par(cex.main=0.8)
heatmap.2( fpkm[ topVarGenes, ], cexCol=0.5, cexRow=0.3, offsetRow=-0.4, offsetCol=-0.4, 
           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), srtCol=30,
           key.par=list(cex=0.6), labRow = lab)
dev.off()


pdf("cummeRbund-heatmap-topVarGenes.pdf", w=8, h=8)
heatmap.2( fpkm[ topVarGenes, ], cexCol=0.9, cexRow=0.5, offsetRow=-0.4, offsetCol=-0.4, 
           scale="row", trace="none", dendrogram="none", main="Top 100 Variance Genes Heatmap",
           Colv=FALSE, col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), labRow = lab, 
           srtCol=30)
dev.off()

# Sample distance heatmap
mycols <- brewer.pal(8, "Dark2")[1:length(unique(DataGroups))]
sampleDists <- as.matrix(dist(t(fpkm)))
png("cummeRbund-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

pdf("cummeRbund-heatmap-samples.pdf", w=8, h=8)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[DataGroups], RowSideColors=mycols[DataGroups],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()
