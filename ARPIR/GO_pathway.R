#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

table_path <- args[1]
dea <- args[2]
output_dir <- args[3]
max_c <- args[4]
script <- args[5]

setwd(output_dir)
max_c = as.numeric(max_c)

#####################################################
######## Gene Ontology enrichment analysis ##########
#####################################################

print("######### Gene Ontology #########")

if (dea == "edgeR"|dea == "DESeq2") {
  table=read.csv(table_path,row.names = 1)
  table1=table[,c(1,2,4)]
}

if (dea == "cummeRbund") {
  table=read.csv(table_path,row.names = 1)
  table1=table[,c(1,2,3)]
}

colnames(table1)=c("Gene","logFC","FDR")
table1=table1[table1$logFC!=Inf&table1$logFC!=-Inf,]

library(rlist)
library(clusterProfiler)
library(org.Hs.eg.db)
library(igraph)
library(scales)
source(paste0(script,"/CNETPLOT_FUNCTION.R"))

# GO tables
table1$entrez = mapIds(org.Hs.eg.db, keys=as.character(table1$Gene), column="ENTREZID",
                       keytype="SYMBOL", multiVals="first")
GO_BP=enrichGO(table1[abs(table1$logFC)>1.5 & table1$FDR<0.05,"entrez"], ont="BP", pvalueCutoff= 0.05, 
               pAdjustMethod = "fdr", universe = table1$entrez, OrgDb = org.Hs.eg.db, readable = T)
GO_CC=enrichGO(table1[abs(table1$logFC)>1.5 & table1$FDR<0.05,"entrez"], ont="CC", pvalueCutoff= 0.05, 
               pAdjustMethod = "fdr", universe = table1$entrez, OrgDb = org.Hs.eg.db, readable = T)
GO_MF=enrichGO(table1[abs(table1$logFC)>1.5 & table1$FDR<0.05,"entrez"], ont="MF", pvalueCutoff= 0.05, 
               pAdjustMethod = "fdr", universe = table1$entrez, OrgDb = org.Hs.eg.db, readable = T)
res_all=list()
if (nrow(GO_BP@result)!=0) {
  GO_BP@result$GO_domain="biological_process"
  res_all=list.append(res_all,GO_BP@result)
}
if (nrow(GO_CC@result)!=0) {
  GO_CC@result$GO_domain="cellular_component"
  res_all=list.append(res_all,GO_CC@result)
}
if (nrow(GO_MF@result)!=0) {
  GO_MF@result$GO_domain="molecular_function"
  res_all=list.append(res_all,GO_MF@result)
}
GO=rbind(GO_BP@result,GO_CC@result,GO_MF@result)
write.csv(GO,paste0(output_dir,"/Gene_ontology/GO_fc1.5_pv0.05.csv"),row.names = F)

# treemap
print("#### Treemap ####")
library(treemap)
tree_all=list()
for (i in 1:length(res_all)) {
  tree=res_all[[i]][,c(2,9)]
  tree$namespace_1003=res_all[[i]][1,10]
  tree_all=list.append(tree_all,tree)
}
all=data.frame(matrix(ncol = ncol(tree_all[[1]]), nrow = 0))
colnames(all) <- colnames(tree_all[[1]])
for (i in 1:length(tree_all)) {all=merge(all,tree_all[[i]],all.x=T,all.y=T)}
pdf(paste0(output_dir,"/Gene_ontology/treemap_GO_fc1.5_pv0.05.pdf"),12,10)
treemap(all,index=c("namespace_1003","Description"),vSize="Count",type="categorical",
        inflate.labels = TRUE,vColor = "namespace_1003",fontsize.labels = c(0,0.9),
        title="GO |Fold-Change|> 1.5 and pvalue<.05",title.legend = "Legend",
        palette=c("lightpink","skyblue","lightgreen"))
dev.off()
png(paste0(output_dir,"/Gene_ontology/treemap_GO_fc1.5_pv0.05.png"),1200, 1000, pointsize=20)
treemap(all,index=c("namespace_1003","Description"),vSize="Count",type="categorical",
        inflate.labels = TRUE,vColor = "namespace_1003",fontsize.labels = c(0,0.9),
        title="GO |Fold-Change|> 1.5 and pvalue<.05",title.legend = "Legend",
        palette=c("lightpink","skyblue","lightgreen"))
dev.off()

# GO tables for genes
print("#### GO tables for genes ####")
GO_gene=list()
for (j in 1:length(res_all)) {
  gene=c()
  go=c()
  a=res_all[[j]]$geneID
  for (i in 1:length(a)) {
    b=unlist(strsplit(a[i], split="/"))
    gene=c(gene,b)
    go=c(go,rep(res_all[[j]]$Description[i],length(b)))
  }
  tab=data.frame(Gene=gene,GO=go)
  y=table1[,c(1,2,3)]
  tab1=merge(tab,y,all.x=T,all.y=F)
  tab1=tab1[order(tab1$GO),]
  tab1$GO_domain=res_all[[j]][1,10]
  GO_gene=list.append(GO_gene,tab1)
}
tab_all=do.call("rbind", GO_gene)
write.csv(tab_all,paste0(output_dir,"/Gene_ontology/tab_GO_genes.csv"),row.names = F)

# dotplot
print("#### Dotplot ####")
if (nrow(GO_BP@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/dotplot_GO_BP.pdf"),10,8)
  print(dotplot(GO_BP, showCategory=max_c, title = "Enriched GO for Biological Process"))
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/dotplot_GO_BP.png"),1000, 800, pointsize=20)
  print(dotplot(GO_BP, showCategory=max_c, title = "Enriched GO for Biological Process"))
  dev.off()
}

if (nrow(GO_CC@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/dotplot_GO_CC.pdf"),10,8)
  print(dotplot(GO_CC, showCategory=max_c, title = "Enriched GO for Cellular Component"))
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/dotplot_GO_CC.png"),1000, 800, pointsize=20)
  print(dotplot(GO_CC, showCategory=max_c, title = "Enriched GO for Cellular Component"))
  dev.off()
}

if (nrow(GO_MF@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/dotplot_GO_MF.pdf"),10,8)
  print(dotplot(GO_MF, showCategory=max_c, title = "Enriched GO for Molecular Function"))
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/dotplot_GO_MF.png"),1000, 800, pointsize=20)
  print(dotplot(GO_MF, showCategory=max_c, title = "Enriched GO for Molecular Function"))
  dev.off()
}

#cnetplot
print("#### Cnetplot ####")
prov=y[order(y$logFC,decreasing = T),]
prov=prov[!duplicated(prov$Gene), ]
prov=prov[!is.na(prov$Gene), ]
rownames(prov)=prov$Gene

geneList=prov$logFC
names(geneList)=prov$Gene

min=-(round(((max(abs(geneList))+1)-1),1))
max=round(((max(abs(geneList))+1)+1),1)

if (nrow(GO_BP@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/cnetplot_GO_BP.pdf"),10,8)
  cnetplot.enrichResult(GO_BP,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Biological Process FC>1.5 pV<.05",
                        cex.main=1)
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/cnetplot_GO_BP.png"),1000, 800, pointsize=20)
  cnetplot.enrichResult(GO_BP,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Biological Process FC>1.5 pV<.05")
  dev.off()
}

if (nrow(GO_CC@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/cnetplot_GO_CC.pdf"),10,8)
  cnetplot.enrichResult(GO_CC,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Cellular Component FC>1.5 pV<.05")
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/cnetplot_GO_CC.png"),1000, 800, pointsize=20)
  cnetplot.enrichResult(GO_CC,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Cellular Component FC>1.5 pV<.05")
  dev.off()
}

if (nrow(GO_MF@result)!=0) {
  pdf(paste0(output_dir,"/Gene_ontology/cnetplot_GO_MF.pdf"),10,8)
  cnetplot.enrichResult(GO_MF,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Molecular Function FC>1.5 pV<.05")
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/cnetplot_GO_MF.png"),1000, 800, pointsize=20)
  cnetplot.enrichResult(GO_MF,categorySize="pvalue", foldChange=geneList, showCategory = max_c,
                        col.bin=seq(min, max, by = 1), main="GO Molecular Function FC>1.5 pV<.05")
  dev.off()
}


##############################################
######## Pathway enrichment analysis #########
##############################################

print("########### Pathway enrichment analysis  ###########")

gene=table1[abs(table1$logFC)>1.5 & table1$FDR<0.05,"entrez"]
tmp=table1[order(table1$logFC,decreasing = T),]
genelist=tmp$logFC
names(genelist)=tmp$Gene

kk=enrichKEGG(gene=gene, organism='hsa', pvalueCutoff = 0.05, pAdjustMethod = "fdr", 
              universe=table1$entrez)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db,keytype = "ENTREZID")
write.csv(kk@result,paste0(output_dir,"/Pathway_analysis/pathway_FC1.5_pv0.05.csv"),row.names = F)

#dotplot
print("#### Dotplot ####")
pdf(paste0(output_dir,"/Pathway_analysis/dotplot_pathways.pdf"),10,8)
dotplot(kk,showCategory=max_c, title = "Enriched pathways")
dev.off()
png(paste0(output_dir,"/Pathway_analysis/dotplot_pathways.png"),1000, 800, pointsize=20)
dotplot(kk,showCategory=max_c, title = "Enriched pathways")
dev.off()

#cnetplot
print("#### Cnetplot ####")
pdf(paste0(output_dir,"/Pathway_analysis/cnetplot_pathways.pdf"),10,8)
par(cex.main=1)
cnetplot.enrichResult(kk, categorySize="pvalue", foldChange=genelist, showCategory = max_c,
                      main="Enriched pathways FC>1.5 & pV<.05", col.bin=seq(min, max, by = 1))
dev.off()
png(paste0(output_dir,"/Pathway_analysis/cnetplot_pathways.png"),1000, 800, pointsize=20)
cnetplot.enrichResult(kk, categorySize="pvalue", foldChange=genelist, showCategory = max_c,
                      main="Enriched pathways FC>1.5 & pV<.05", col.bin=seq(min, max, by = 1))
dev.off()

#pathview
print("#### Pathview ####")
library(pathview)
new_data=table1[,c("logFC","entrez")]
new_data=new_data[!duplicated(new_data$entrez), ]
new_data=new_data[!is.na(new_data$entrez), ]
rownames(new_data)=new_data$entrez
new_data=new_data[,-2,drop=F]
setwd(paste0(output_dir,"/Pathway_analysis/pathview"))
tmp = sapply(kk@result$ID, function(pid) tryCatch(pathview(gene.data=new_data, pathway.id=pid, species="hsa",
                                                           low="dodgerblue",high="firebrick1",mid="gray88"),error=function(e) NULL))


# histogram go
print("#### Histogram ####")
library(dplyr)
library(ggplot2)

for (j in 1:length(unique(tab_all$GO_domain))) {
  if (nrow(GO[GO$GO_domain==unique(tab_all$GO_domain)[j],])>=30) {
    names=as.character(GO$Description[GO$GO_domain==unique(tab_all$GO_domain)[j]][1:30])
  } else {
    names=as.character(GO$Description[GO$GO_domain==unique(tab_all$GO_domain)[j]][1:nrow(GO[GO$GO_domain==unique(tab_all$GO_domain)[j],])])
  }
  p=tab_all[tab_all$GO %in% names,]
  p1=p[p$FDR<0.05 & abs(p$logFC)>1.5,]
  p1$fc=ifelse(p1$logFC>=0,"up","down")
  p1=p1[!is.na(p1$Gene),c("Gene","GO","fc","logFC")]
  
  plotting_df <-
    p1 %>% 
    group_by(GO, fc) %>% 
    summarise(Freq = n()) %>% 
    mutate(Freq = if_else(fc == "down", -Freq, Freq))
  plotting_df$FoldChange=NA
  for (i in 1:nrow(plotting_df)) {
    plotting_df$FoldChange[i]=mean(p1$logFC[p1$GO==plotting_df$GO[i] & 
                                              p1$fc==plotting_df$fc[i]])
  }
  temp_df <-
    plotting_df %>% 
    filter(fc == "up") %>% 
    arrange(Freq)
  the_order <- temp_df$GO
  q <- 
    plotting_df %>% 
    ggplot(aes(x = GO, y = Freq, fill= FoldChange)) +
    geom_bar(stat = "identity", width = 0.75) +
    coord_flip() +
    scale_x_discrete(limits = the_order) +
    scale_y_continuous(breaks = seq(-300, 300, 10), labels = abs(seq(-300, 300, 10))) +
    labs(x = "GO", y = "Gene count", title = "\nGene Ontology", fill="FoldChange mean",
         subtitle = paste0(unique(tab_all$GO_domain)[j],"\n")) +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5,face = "bold", size = 16),
          plot.subtitle = element_text(hjust = 0.5,face = "bold", size = 12),
          panel.background = element_rect(fill =  "grey90")) +
    scale_fill_gradient2(midpoint=0, low="dodgerblue1", mid="white", high="firebrick2")
  
  pdf(paste0(output_dir,"/Gene_ontology/hist_",unique(tab_all$GO_domain)[j],".pdf"),10,6)
  print(q)
  dev.off()
  png(paste0(output_dir,"/Gene_ontology/hist_",unique(tab_all$GO_domain)[j],".png"),width = 10, height = 6, units = 'in', res = 300)
  print(q)
  dev.off()
}

# histogram path
if (nrow(kk@result)>=30) {
  names=as.character(kk@result$Description[1:30])
} else {
  names=as.character(kk@result$Description[1:nrow(kk@result)])
}
p=list()
for (i in 1:length(names)) {
  p=list.append(p,table1[table1$Gene %in% 
                           strsplit(as.character(
                             kk@result$geneID[kk@result$Description==names[i]]),"/")[[1]],])
}
for (i in 1:length(p)) {p[[i]]$path=names[i]}
p1=data.frame()
for (i in 1:length(p)) {p1=rbind(p1,p[[i]])}
p_go=p1[p1$FDR<0.05 & abs(p1$logFC)>1.5,]
p_go$fc=ifelse(p_go$logFC>=0,"up","down")
p_go=p_go[!is.na(p_go$Gene),c("Gene","path","fc","logFC")]
plotting_df <-
  p_go %>% 
  group_by(path, fc) %>% 
  summarise(Freq = n()) %>% 
  mutate(Freq = if_else(fc == "down", -Freq, Freq))
plotting_df$FoldChange=NA
for (i in 1:nrow(plotting_df)) {
  plotting_df$FoldChange[i]=mean(p_go$logFC[p_go$path==plotting_df$path[i] & 
                                              p_go$fc==plotting_df$fc[i]])
}
temp_df <-
  plotting_df %>% 
  filter(fc == "up") %>% 
  arrange(Freq)
the_order <- temp_df$path
q <- 
  plotting_df %>% 
  ggplot(aes(x = path, y = Freq, fill= FoldChange)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = the_order) +
  scale_y_continuous(breaks = seq(-300, 300, 10), 
                     labels = abs(seq(-300, 300, 10))) +
  labs(x = "Pathway", y = "Gene count", title = "\nPathways\n", fill="FoldChange mean") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5,face = "bold", size = 16),
        panel.background = element_rect(fill =  "grey90")) +
  scale_fill_gradient2(midpoint=0, low="dodgerblue1", mid="white", high="firebrick2")

pdf(paste0(output_dir,"/Pathway_analysis/hist_pathway.pdf"),10,6)
print(q)
dev.off()
png(paste0(output_dir,"/Pathway_analysis/hist_pathway.png"),width = 10, height = 6, units = 'in', res = 300)
print(q)
dev.off()
