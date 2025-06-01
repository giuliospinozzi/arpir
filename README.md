# ARPIR (Automatic RNA-Seq Pipelines with Interactive Report)

## Introduction
<p align="justify"> ARPIR makes RNA-Seq analysis: quality control, pre-processing, alignment (primary-analysis), transcript quantification and differential expression analysis on BAM files (secondaty-analysis). Given the input files and the working directory, the pipeline is completely automated. First, quality control on FastQ files is performed with FastQC and FastQ Screen. FastQC makes quality control and creates one report for sample. FastQ Screen estimates approximately the percentage of reads that can be mapped on genomes other than human, like ribosomal genomes, PhiX genome and mouse genome. This allows to evaluate the presence of contaminating genomes. Pre-processing follows quality control: the reads are aligned on PhiX genome and ribosomal genome to eliminate contaminations. Alignment can be performed with TopHat2, HISAT2 or STAR; in the first case quantification is performed with Cufflinks and DEA with cummeRbund, in the second case quantification is performed with featureCounts and DEA with DESeq2 or edgeR, in the third case it is possible to choose one of the previous methods for quantification and DEA. A second intermediate quality control analysis is also performed on the aligned BAM files with some of the RSeQC scripts and in particular: inner_distance, junction_annotation, junction_saturation, bam_stat, read_distribution, geneBody_coverage. A final report is generated from the results using multiqc. It is possible to perform an optional tertiary-analysis on the results. It consists in Gene Ontology enrichment analysis and KEGG Pathway enrichment analysis on the differentially expressed genes (with absolute Fold Change value higher than 1.5 and adjusted p-value lower than 0.05). </p>

## Install Option 1 (suggested)
If you have a newer system please use Conda to install some important packages and comment in the configuration script the relative lines:

```
conda env create -f arpir.yml

```
Change py interpreter in each py script used (/opt/anaconda3/envs/arpir/bin/python2)

## Install Option 2 (old)
<p align="justify"> For ARPIR to work properly, you must first make sure that you have installed the necessary applications for the pipeline. You can use the conf.sh configuration file to install all the required software and dependencies, as well as download the genome files. The script runs on Debian-based systems and you can run it with the following command: </br>

```sudo bash conf.sh```

## Install Option 3 (manually)

If instead you want to manually proceed with the installation of the prerequisites, the necessary ones are the following: </p>

### Applications
#### For Graphical User Interface and scripts running:
*	zenity 3.18.1.1 (https://help.gnome.org/users/zenity/)
*	Python 2.7.12 (modules: os, argparse, sys, csv, pandas) (https://www.python.org/downloads/)
*	R 3.4.3 (packages: cummeRbund, edgeR, DESeq2, ggfortify, ggrepel, genefilter, RColorBrewer, gplots, clusterProfiler, dplyr, org.Hs.eg.db, igraph, scales, treemap, pathview, shiny, DT, magick, rlist, visNetwork, shinyjs, knitr) (https://www.r-project.org/)
* pandoc (https://pandoc.org/)
* multiqc 1.7 (https://multiqc.info/)
#### For quality control on fastq files and pre-processing:
*	FastQC 0.11.5 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
*	FastQ Screen 0.11.3 (https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
*	bwa 0.7.12-r1039 (http://bio-bwa.sourceforge.net/)
*	samtools 1.9 (http://samtools.sourceforge.net/)
* pigz
* GD::Graph Perl module
#### For alignment and quality control on BAM files:
*	TopHat 2.1.1 (https://ccb.jhu.edu/software/tophat/index.shtml)
*	RSeQC 2.6.4 (http://rseqc.sourceforge.net/)
*	HISAT 2.1.0 (https://ccb.jhu.edu/software/hisat2/index.shtml)
* STAR 2.5.0a (https://github.com/alexdobin/STAR)
#### For quantification:
*	featureCounts 1.5.3 (http://subread.sourceforge.net/)
*	Cufflinks 2.2.1 (http://cole-trapnell-lab.github.io/cufflinks/)
### Input files
The necessary input files are as follows:
*	<p align="justify">FastQ files of the samples. Since this is a differential analysis, at least two cases and two controls must be present, furthermore, if it is a paired-end analysis, files of read 1 and read 2 must be provided for each sample while if it is a single-end analysis only files of read 1 must be provided.</p>
*	<p align="justify">Reference genome for Bowtie2 with relative indexes, necessary for alignment with TopHat2. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/index/bowtie2/hg19. It is not necessary to report the extension, but only the name of the file and the indexes must be present in the same folder.</p>
*	<p align="justify">Reference genome with relative indexes, necessary for alignment with HISAT2. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/index/hisat2/hg19. The reference genome file is the same as TopHat, but the indexes are obtained differently (see below).</p>
* <p align="justufy">Reference genome with relative indexes, necessary for alignment with STAR. The path in which the directory with indexes is searched by default is as follows: /opt/genome/human/hg19/index/STAR. The reference genome file is the same as TopHat2 and HISAT2, but the indexes are obtained differently (see below).</p>
*	<p align="justify">BED file, necessary for quality control on BAM files. The path in which the file is searched by default is as follows: /opt/genome/human/hg19/annotation/hg19.refseq.bed12.</p>
*	<p align="justify">PhiX genome with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/control/phix174/bwa/phiX174.fa.</p>
*	<p align="justify">Ribosomal genome 5s with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa.</p>
*	<p align="justify">Ribosomal genome with relative indexes, necessary for pre-processing. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/contam/bwa/humRibosomal.fa.</p>
*	<p align="justify">GTF file, necessary for alignment with TopHat2 and quantification with featureCounts and Cufflinks. The path in which the file is searched by default is as follows: /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf.</p>
*	<p align="justify">Reference genome with relative indexes, necessary for quantification with Cufflinks. The path in which the file with its indexes is searched by default is as follows: /opt/genome/human/hg19/index/hg19.fa. The reference genome file is the same, but the indexes change. </p>

#### Creating indexes:
<p align="justify"> The indexes can be created starting from the Fasta file and require different commands based on the alignment software for which they will be used.</p>
In the case of TopHat, the necessary command is the following:<br>

```bowtie2-build file.fa file``` <br>

In the case of HISAT, the necessary command is:<br>

```hisat2-build file.fa file```

In the case of STAR, the necessary command is:

```STAR --runMode genomeGenerate --genomeDir index_directory --genomeFastaFiles file.fa --sjdbGTFfile file.gtf```

In the case of BWA, instead, the necessary command is:<br>

```bwa index file.fa```

Finally, in the case of BWA, the necessary command is:<br>

```samtools faidx file.fa```

## Running application
### Command line
<p align="justify"> After installing the necessary tools and getting the required input files, you can launch the application using command-line. The following is the use for analysis in single-end and paired-end, with the arguments that must necessarily be inserted, as they do not provide default options, and the optional arguments, which can be omitted if the paths and options you want match the standard ones.</p>

#### Usage

#### Paired-end

```python ARPIR.py -n <project_name> -pn <pool_name> -sn <sample_A,sample_B,sample_C,sample_D> -r1 <sample_A_read1,sample_B_read1,sample_C_read1,sample_D_read1> -r2 <sample_A_read2,sample_B_read2,sample_C_read2,sample_D_read2> -type <cntrl,cntrl,treat,treat> -comp <cntrl_VS_treat> -o <output_directory> [options]*```

#### Single-end

```python ARPIR.py -n <project_name> -pn <pool_name> -sn <sample_A,sample_B,sample_C,sample_D> -r1 <sample_A_read1,sample_B_read1,sample_C_read1,sample_D_read1> -type <cntrl,cntrl,treat,treat> -comp <cntrl_VS_treat> -o <output_directory> [options]*```

#### Arguments
| | |
------------ | -------------
```-n```	| Project name. No default options. <br>
```-pn```	| Pool name. No default options. <br>
```-sn```	| Sample names (',' sep). No default options. <br>
```-r1```	| Read 1 FastQ path (',' sep). No default options. Files must appear in the same order as sample names. <br>
```-r2```	| Read 2 FastQ path (',' sep). Required only for paired-end analysis. Files must appear in the same order as sample names. <br>
```-type```	| Sample types (',' sep). No default options. Types must appear in the same order as sample names. <br>
```-comp```	| Comparisons between samples (allows multiple comparisons comma separated). No default options. <br>
```-o``` | Output directory. No default options. <br>

#### Options
| | |
------------ | -------------
```-rb```	| Reference genome file path for bowtie. Default: `/opt/genome/human/hg19/index/bowtie2/hg19` <br>
```-rh```	| Reference genome file path for HISAT2. Default: `/opt/genome/human/hg19/index/hisat2/hg19` <br>
```-rs```	| Index directory path for STAR. Default: `/opt/genome/human/hg19/index/STAR` <br>
```-bed```	| Reference genome annotation file path. Default: `/opt/genome/human/hg19/annotation/hg19.refseq.bed12` <br>
```-ph```	| PhiX genome file path. Default: `/opt/genome/control/phix174/bwa/phiX174.fa` <br>
```-rib1```	| Ribosomal genome file path. Default: `/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa` <br>
```-rib2```	| Ribosomal genome file path. Default: `/opt/genome/human/hg19/contam/bwa/humRibosomal.fa` <br>
```-t```	| Max thread number. Default: 12. <br>
```-g```	| GTF file path. Default: `/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf` <br>
```-a```	| Alignment method. Default: hisat; alternative: tophat, star. <br>
```-l```	| Library type. Default: fr-firststrand; alternative: fr-secondstrand. <br>
```-q```	| Quantification method. Default: featureCounts; alternative: Cufflinks. <br>
```-r```	| Reference genome file path. Default: `/opt/genome/human/hg19/index/hg19.fa` <br>
```-dea```	| Differential Expression Analysis method. Default: edgeR; alternatives: DESeq2, cummeRbund. <br>
```-r_path```	| Script directory (alignment, quantification and DEA). Default: `/opt/applications/src/arpir/ARPIR`. <br>
```-meta``` | Analysis with or without final tertiary-analysis. Default: full; alternative: quant. <br>
```-cat``` | Max number of categories showed in R plots for tertiary-analysis. Default: 5. <br>

#### Example
<p align="justify"> To make the launch of the command easier, we have created a script run.sh where you can enter all the desired variables. Once done, you can launch the pipeline with the following command:</p>

```bash run.sh```

### Graphical User Interface
#### Usage
```bash GUI_ARPIR.sh``` <br>
<p align="justify"> A series of windows allows to indicate input files and to choose between various options.</p>

#### Example
```bash /path/to/script/ARPIR/GUI_ARPIR.sh ``` <br>
<p align="justify"> A first window will appear in which you will be asked to enter the data relating to the experiment, as in the example below.</p>

<p align="center"><img src="/images/GUI_screenshot1.png" width="55%"></p>

<p align="justify"> This is followed by a series of windows that allow you to choose between different options or ask to select the input files or folders, as shown in the following two images.</p>

<p align="center"><img src="/images/GUI_screenshot2.png" width="45%"></p>

<p align="center"><img src="/images/GUI_screenshot3.png" width="25%"></p>

<p align="justify"> If they exist, the default paths will be shown automatically. Finally, once all the inputs and the parameters have been selected, a summary window will appear. If the data entered are correct, it will be sufficient to confirm to start the analysis, otherwise it is possible to cancel and start from the beginning.</p>

### Shiny app
To observe the results of the analysis in an interactive form, it is possible to launch a [shiny app](shiny.md).

