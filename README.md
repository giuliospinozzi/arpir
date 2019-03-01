# ARPIR (Automatic RNA-Seq Pipelines with Interactive Report)

## Introduction
<p align="justify"> This application makes RNA-Seq analysis: quality control, pre-processing, alignment, transcript quantification and differential expression analysis on BAM files. Given the input files and the working directory, the pipeline is completely automated. First, quality control on fastq files is performed with FastQC e FastQ Screen. FastQC makes quality control and creates one report for sample. FastQ Screen estimates approximately the percentage of reads that can be mapped on genomes other than human, like ribosomal genome, phix genome and mouse genome. This allows to evaluate the presence of contaminating genomes. Pre-processing follows quality control: the reads are aligned on phix genome and ribosomal genome to eliminate contaminations. Alignment can be performed with TopHat or HISAT2; in the first case quantification is performed with Cufflinks and DEA with cummeRbund, in the second case quantification is performed with featureCounts and DEA with DESeq2 or edgeR. A second intermediate quality control analysis is also performed on the aligned BAM files with some of the RSeQC scripts and in particular: inner_distance, junction_annotation, junction_saturation, bam_stat, read_distribution, geneBody_coverage. It is possible to perform an optional meta-analysis on the results. It consists in Gene Ontology enrichment analysis and KEGG Pathway enrichment analysis on the differentially expressed genes (with absolute Fold Change value higher than 1.5 and adjusted p-value lower than 0.05). </p>

## Prerequisites
### Applications
#### For Graphical User Interface and scripts running:
*	zenity 3.18.1.1
*	Python 2.7.12 (modules: os, argparse, sys, csv, pandas)
*	R 3.4.3 (packages: cummeRbund, edgeR, DESeq2, ggfortify, ggrepel, genefilter, RColorBrewer, gplots, clusterProfiler, dplyr, org.Hs.eg.db, igraph, scales, treemap, pathview, shiny, DT, magick, rlist, visNetwork, shinyjs, knitr)
* pandoc
#### For quality control on fastq files and pre-processing:
*	FastQC 0.11.5
*	FastQ Screen 0.11.3
*	bwa 0.7.12-r1039
*	samtools 0.1.19-96b5f2294a
*	mysql 14.14
#### For alignment and quality control on BAM files:
*	TopHat 2.1.1
*	RSeQC 2.6.4
*	HISAT 2.1.0
#### For quantification:
*	featureCounts 1.5.3
*	Cufflinks 2.2.1
### Input files
*	<p align="justify">Fastq files of the samples. Since this is a differential analysis, at least two case and two control must be present, furthermore if it is a paired-end analysis files of read 1 and read 2 must be provided for each sample while if it is a single-end analysis only files of read 1 must be provided.</p>
*	<p align="justify">Reference genome for Bowtie with relative indexes, necessary for alignment with TopHat. The path in which the file with its indices is searched by default is as follows: /opt/genome/human/hg19/index/bowtie2/hg19.</p>
*	<p align="justify">Reference genome for Hisat2 with relative indexes, necessary for alignment with HISAT2. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/index/hisat2/hg19.</p>
*	<p align="justify">BED file, necessary for quality control on BAM files. The path in which the file is searched by default is as follows: /opt/genome/human/hg19/annotation/hg19.refseq.bed12.</p>
*	<p align="justify">Phix genome with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/control/phix174/bwa/phiX174.fa.</p>
*	<p align="justify">Ribosomal genome 1 with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa.</p>
*	<p align="justify">Ribosomal genome 2 with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/contam/bwa/humRibosomal.fa.</p>
*	<p align="justify">GTF file, necessary for alignment with TopHat and quantification with featureCounts and Cufflinks. The path in which the file is searched by default is as follows: /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf. </p>
*	<p align="justify">Reference genome with relative indexes, necessary for quantification with Cufflinks. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/index/hg19.fa. </p>

## Running application
### Command line
#### Usage

#### Paired-end

```python ARPIR.py -n <project_name> -pn <pool_name> -sn <sample_A,sample_B,sample_C,sample_D> -r1 <sample_A_read1,sample_B_read1,sample_C_read1,sample_D_read1> -r2 <sample_A_read2,sample_B_read2,sample_C_read2,sample_D_read2> -type <cntrl,cntrl,treat,treat> -o <output_directory> [options]*```

#### Single-end

```python ARPIR.py -n <project_name> -pn <pool_name> -sn <sample_A,sample_B,sample_C,sample_D> -r1 <sample_A_read1,sample_B_read1,sample_C_read1,sample_D_read1> -type <cntrl,cntrl,treat,treat> -o <output_directory> [options]*```

#### Arguments
| | |
------------ | -------------
```-n```	| Project name. No default option. <br>
```-pn```	| Pool name. No default option. <br>
```-sn```	| Sample names (',' sep). No default option. <br>
```-r1```	| Read 1 fastq path (',' sep). No default option. Files must appear in the same order as sample names. <br>
```-r2```	| Read 2 fastq path (',' sep). Required only for paired-end analysis. Files must appear in the same order as sample names. <br>
```-type```	| Sample types (',' sep, 'cntrl' for control). No default option. Types must appear in the same order as sample names. <br>
```-o``` | Output directory. No default option. <br>

#### Options
| | |
------------ | -------------
```-rb```	| Reference genome file path for bowtie. Default: `/opt/genome/human/hg19/index/bowtie2/hg19` <br>
```-rh```	| Reference genome file path for hisat2. Default: `/opt/genome/human/hg19/index/hisat2/hg19` <br>
```-bed```	| Reference genome annotation file path. Default: `/opt/genome/human/hg19/annotation/hg19.refseq.bed12` <br>
```-ph```	| Phix genome file path. Default: `/opt/genome/control/phix174/bwa/phiX174.fa` <br>
```-rib1```	| Ribosomal genome file path. Default: `/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa` <br>
```-rib2```	| Ribosomal genome file path. Default: `/opt/genome/human/hg19/contam/bwa/humRibosomal.fa` <br>
```-t```	| Max thread number. Default: 12. <br>
```-g```	| GTF file path. Default: `/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf` <br>
```-a```	| Alignment method. Default: hisat; alternative: tophat. <br>
```-l```	| Library type. Default: fr-firststrand; alternative: fr-secondstrand. <br>
```-q```	| Quantification method. Default: featureCounts; alternative: Cufflinks. <br>
```-r```	| Reference genome file path. Default: `/opt/genome/human/hg19/index/hg19.fa` <br>
```-dea```	| Differential Expression Analysis method. Default: edgeR; alternatives: DESeq2, cummeRbund. <br>
```-meta``` | Analysis with or without final meta-analysis. Default: full; alternative: quant. <br>
```-cat``` | Max number of category showed in R plots for meta-analysis. Default: 5. <br>

### Graphical User Interface
#### Usage
```bash GUI_ARPIR.sh``` <br>
A series of windows allow to indicate input files and to choose between various options.

<p align="center"><img src="/images/GUI_screenshot.png" width="55%"></p>

### Shiny app
To observe the results of the analysis in an interactive form, it is possible to launch a [shiny app](shiny.md).

