######## BALLGOWN #########
library(ballgown)

pheno_data=data.frame(ids=c("sample_01","sample_02","sample_03","sample_04","sample_05","sample_06"),
                      type=c("Cntrl","Cntrl","Cntrl","Treat","Treat","Treat"),
                      path=c("Hisat/StringTie/sample_01","Hisat/StringTie/sample_02",
                             "Hisat/StringTie/sample_03","Hisat/StringTie/sample_04",
                             "Hisat/StringTie/sample_05","Hisat/StringTie/sample_06"))

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples=as.vector(pheno_data$path), pData=pheno_data)

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))
results_genes=results_genes[results_genes$pval!="NaN",]
write.csv(results_genes,"ballgown.csv")

######## CUMMERBUND #########
library(cummeRbund)
cuff_data <- readCufflinks("./TopHat/cuffdiff/",rebuild = T)

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
diffGenesOutput1=diffGenesOutput[,c(1,7,9)]

diffGenesOutput2=diffGenesOutput1[diffGenesOutput1$p_value<1,]

for (i in 1:nrow(diffGenesOutput2)) {
  diffGenesOutput2$gene[i]=strsplit(as.character(diffGenesOutput2$x[i]),",")[[1]][1]
}
write.csv(diffGenesOutput2,"cummerbund.csv")


######## DESEQ2 #########
# Import data from featureCounts
countdata1 <- read.table("./Hisat/featureCounts_counts.txt",header=TRUE, row.names=1,sep = "\t")
lenth_genes <- countdata1[ ,5,drop=F]
countdata <- countdata1[ ,6:ncol(countdata1)]
colnames(countdata)=c("sample_01","sample_02","sample_03","sample_04","sample_05","sample_06")
DataGroups <- factor(c("cntrl","cntrl","cntrl","treat","treat","treat"))
DataGroups <- relevel(DataGroups,ref="cntrl")

# Convert to matrix
countdata <- as.matrix(countdata)

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet
coldata <- data.frame(row.names=colnames(countdata), DataGroups)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~DataGroups)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Get differential expression results
res <- results(dds)
res <- res[order(res$padj), ]
res=as.data.frame(res)
res=res[!is.na(res$padj),]
write.csv(res,"deseq2.csv")


######## EDGER #########
library(edgeR)
dgList <- DGEList(counts=countdata,group=factor(DataGroups),genes = lenth_genes)

# filter data to retain genes that are represented at least 1 counts per million (cpm) in at least 2 samples
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
dgList$samples$lib.size <- colSums(dgList$counts)

# normalization using TMM method
dgList <- calcNormFactors(dgList, method="TMM")

# Dispersion estimates
design.mat <- model.matrix(~ 0 + dgList$samples$group)
colnames(design.mat) <- levels(dgList$samples$group)
dgList <- estimateGLMCommonDisp(dgList,design.mat)
dgList <- estimateGLMTrendedDisp(dgList,design.mat, method="power")
dgList <- estimateGLMTagwiseDisp(dgList,design.mat)

# Differential expression analysis
fit <- glmFit(dgList, design.mat)
lrt <- glmLRT(fit, contrast=c(-1,1))
edgeR_results <- topTags(lrt, n=Inf)
edgeR_results=as.data.frame(edgeR_results)
write.csv(edgeR_results,"edger.csv")


######## SLEUTH #########
library("sleuth")

# from ensembl to gene name
tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}
t2g <- tx2gene()

base_dir <- "./kallisto/"
samples=c("sample_01","sample_02","sample_03","sample_04","sample_05","sample_06")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

s2c <- data.frame(path=kal_dirs, sample=samples, type = rep(c("cntrl","treat"), each=3), stringsAsFactors=FALSE)

# sleuth analysis
so <- sleuth_prep(s2c, ~type, target_mapping = t2g)

# fit the model and compare the two types
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="typetreat") 

save(so,file="sleuth.RData")

# # display the results
# sleuth_live(so)
# 
# res=read.csv("./kallisto/test_table.csv",row.names = 1)
# res=res[!is.na(res$pval),]
# write.csv(res,"kallisto.csv")


##### STAR #####

#### EDGER ####


#### DESEQ2 ####


#### CUMMERBUND ####

