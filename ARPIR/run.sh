# Replace the variables with paths and options of interest.

# Project name.
PROJECT_NAME="Project_01";

# Pool name.
POOL_NAME="Simulation";

# Sample names (',' sep).
SAMPLE_NAMES="sample_01,sample_02,sample_03,sample_04,sample_05,sample_06";

# Read 1 fastq path (',' sep). Files must appear in the same order as sample names.
READ_1="/opt/ngs/Simulation/sample_01_1.fastq.gz,/opt/ngs/Simulation/sample_02_1.fastq.gz,/opt/ngs/Simulation/sample_03_1.fastq.gz,/opt/ngs/Simulation/sample_04_1.fastq.gz,/opt/ngs/Simulation/sample_05_1.fastq.gz,/opt/ngs/Simulation/sample_06_1.fastq.gz";

# Read 2 fastq path (',' sep). Files must appear in the same order as sample names. Required only for paired-end analysis.
READ_2="/opt/ngs/Simulation/sample_01_2.fastq.gz,/opt/ngs/Simulation/sample_02_2.fastq.gz,/opt/ngs/Simulation/sample_03_2.fastq.gz,/opt/ngs/Simulation/sample_04_2.fastq.gz,/opt/ngs/Simulation/sample_05_2.fastq.gz,/opt/ngs/Simulation/sample_06_2.fastq.gz";

# Sample types (',' sep). Types must appear in the same order as sample names.
TYPE="cntrl,cntrl,cntrl,treat,treat,treat";

# Output directory.
OUTPUT_DIR="/opt/ngs";

# Reference genome file path for bowtie.
BOWTIE_DNA="/opt/genome/human/hg19/index/bowtie2/hg19";

# Reference genome file path for hisat2.
HISAT_DNA="/opt/genome/human/hg19/index/hisat2/hg19";

# Reference genome file path for star.
STAR_DNA="/opt/genome/human/hg19/index/STAR";

# Reference genome annotation file path.
BED_FILE="/opt/genome/human/hg19/annotation/hg19.refseq.bed12";

# Phix genome file path.
PHIX_DNA="/opt/genome/control/phix174/bwa/phiX174.fa";

# Ribosomal genome 1 (5S) file path.
RIBOSOMAL_DNA_1="/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa";

# Ribosomal genome 2 file path.
RIBOSOMAL_DNA_2="/opt/genome/human/hg19/contam/bwa/humRibosomal.fa";

# Max thread number.
THREADS="12";

# GTF file path.
GFT_FILE="/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf";

# Alignment method. Options: star, hisat, tophat.
ALIGNMENT="hisat";

# Library type. Options: fr-firststrand, fr-secondstrand.
LIBRARY="fr-firststrand";

# Quantification method. Options: featureCounts, Cufflinks.
QUANTIFICATION="featureCounts";

# Reference genome file path.
REFERENCE_DNA="/opt/genome/human/hg19/index/hg19.fa";

# Differential Expression Analysis method. Options: edgeR, DESeq2, cummeRbund.
DEA="edgeR";

# Script directory.
SCRIPT_DIR="/opt/applications/src/arpir/ARPIR";

# Analysis with or without final meta-analysis. Options: full, quant.
META_ANALYSIS="full";

# Max number of categories showed in R plots for meta-analysis.
CATEGORY_NUMBER="5";

# LOG file.
LOG="/opt/ngs/logs/simulation.log";

#Comparisons (cntrl_VS_treat1,cntrl_VS_treat2).
COMPARISONS="cntrl_VS_treat";


## uncomment the analysis type of interest
# paired-end analysis
#nohup python ${SCRIPT_DIR}/ARPIR.py -n ${PROJECT_NAME} -pn ${POOL_NAME} -sn ${SAMPLE_NAMES} -r1 ${READ_1} -r2 ${READ_2} -type ${TYPE} -o ${OUTPUT_DIR} -rb ${BOWTIE_DNA} -rh ${HISAT_DNA} -rs ${STAR_DNA} -bed ${BED_FILE} -ph ${PHIX_DNA} -rib1 ${RIBOSOMAL_DNA_1} -rib2 ${RIBOSOMAL_DNA_2} -t ${THREADS} -g ${GFT_FILE} -a ${ALIGNMENT} -l ${LIBRARY} -q ${QUANTIFICATION} -r_path ${SCRIPT_DIR} -r ${REFERENCE_DNA} -dea ${DEA} -comp ${COMPARISONS} -meta ${META_ANALYSIS} -cat ${CATEGORY_NUMBER} 2>&1 > ${LOG} &

# single-end analysis
#nohup python ${SCRIPT_DIR}/ARPIR.py -n ${PROJECT_NAME} -pn ${POOL_NAME} -sn ${SAMPLE_NAMES} -r1 ${READ_1} -type ${TYPE} -o ${OUTPUT_DIR} -rb ${BOWTIE_DNA} -rh ${HISAT_DNA} -rs ${STAR_DNA} -bed ${BED_FILE} -ph ${PHIX_DNA} -rib1 ${RIBOSOMAL_DNA_1} -rib2 ${RIBOSOMAL_DNA_2} -t ${THREADS} -g ${GFT_FILE} -a ${ALIGNMENT} -l ${LIBRARY} -q ${QUANTIFICATION} -r_path ${SCRIPT_DIR} -r ${REFERENCE_DNA} -dea ${DEA} -comp ${COMPARISONS} -meta ${META_ANALYSIS} -cat ${CATEGORY_NUMBER} 2>&1 > ${LOG} &
