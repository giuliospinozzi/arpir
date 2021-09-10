#!/bin/bash
source /etc/environment
source /etc/profile

PIPEVERSION="1.0 - creo_RNAseq"
STARTTIME=`date +'%y-%m-%d %H:%M:%S'`
RUN_ID=`date +"%Y%m%d%H%M%S"`

RED='\033[0;31m' 
NC='\033[0m' # No Color
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

echo "
  +--------------------------------------------------------+
  |                                                        |
  |            Illumina Pipeline: CREO RNAseq              |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi, PhD                        |
  |            Valentina Tini                              |
  |  Date:     June 2018                                   |
  |  Contact:  giulio.spinozzi@unipg.it                    |
  |            valy.tini@hotmail.it                        |
  |  Version:  2.0 - CREO - RNAseq alignment               |
  |                  No SAM, HiSeq optimized, Paired-End   |
  |                  TopHat2/HISAT2 - Bowtie2              |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1. PROJECT_NAME [PILOT_STUDY]
	2. POOL_NAME [21092017]
	3. LIBRARY_NAME [library290817A1]
	4. RESULTS_DIR [/opt/ngs/results]
	5. READ1 [R1.fastq.gz]
	6. READ2 [R2.fastq.gz]
	7. MAXTHREADS [12]
	8. REFERENCE_GENOME_BOWTIE [/opt/genome/human/hg19/index/bowtie2/hg19]
	9. REFERENCE_GENOME_HISAT2 [/opt/genome/human/hg19/index/hisat2/hg19]
	10. REFERENCE_GENOME_STAR [/opt/genome/human/hg19/index/STAR]
	11. BED_FILE [/opt/genome/human/hg19/annotation/hg19.refseq.bed12]
	12. PHIX_GENOME [/opt/genome/control/phix174/bwa/phiX174.fa]
	13. RIBOSOMAL_GENOME_hum5SrDNA [/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa]
	14. RIBOSOMAL_GENOME_humRibosomal [/opt/genome/human/hg19/contam/bwa/humRibosomal.fa]
	15. ANALYSIS_PROTOCOL [tophat,hisat,star]
	16. GTF_FILE [/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf]
	17. LIBRARY_TYPE [fr-firststrand,fr-secondstrand]
	18. SAMPLE_TYPE [cntrl,treat]
	19. SCRIPT_DIR [/opt/applications/src/arpir/ARPIR]
"
printf "${GREEN}[CREO] RNAseq Pipeline ${NC}\n"
printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}[CREO] Preprocessing input variables (delimiters:<>)${NC}\n"
## print input variables (check for log utils)
INPUTVARNUM=0
for INPUTVAR in "$@"; do
	let INPUTVARNUM++; 
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
    echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
done

printf "${GREEN}@@@@ Variable Adjustments ${NC}\n"

PROJECT_NAME="${1}";
POOL_NAME="${2}";
LIBRARY_NAME="${3}";
#RAW_DATA_DIR="/opt/ngs/raw_data";
RESULTS_DIR="${4}";
R1_FASTQ="${5}";
R2_FASTQ="${6}";
MAXTHREADS="${7}";
REFERENCE_GENOME_BOWTIE="${8}";
REFERENCE_GENOME_HISAT2="${9}";
REFERENCE_GENOME_STAR="${10}";
BED_FILE="${11}";
PHIX_GENOME="${12}";
RIBOSOMAL_GENOME_hum5SrDNA="${13}";
RIBOSOMAL_GENOME_humRibosomal="${14}";
ANALYSIS_PROTOCOL="${15}";
GTF_FILE="${16}";
LIBRARY_TYPE="${17}";
SAMPLE_TYPE="${18}";
SCRIPT_DIR="${19}";

printf "${GREEN}@@@@ Folder creation --> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}${NC}\n"

mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/Quality

RUN_NAME="${PROJECT_NAME}|${POOL_NAME}|${LIBRARY_NAME}"


NUMBER_RAW_READS=$((`zcat ${R1_FASTQ} | wc -l`/4)) ;
printf "${GREEN}@@@@ Starting Number of Raw Reads --> ${NUMBER_RAW_READS}${NC}\n"

# printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Converting BCL in FASTQ files (compressed) #####${NC}\n"
# bcl2fastq --runfolder-dir ${RAW_DATA_DIR}/${PROJECT_NAME}/170913_D00793_0018_Ahv7khbcxy -p 12 --output-dir /opt/ngs/raw_data/${PROJECT_NAME}/170913_D00793_0018_Ahv7khbcxy/FASTQ

printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### FastQC Report #######${NC}\n"
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/Quality/${LIBRARY_NAME}
fastqc --nogroup --extract -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/Quality/${LIBRARY_NAME} -t ${MAXTHREADS} -f fastq ${R1_FASTQ} ${R2_FASTQ}

printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### FastQ Screen Report #######${NC}\n"
fastq_screen ${R1_FASTQ} ${R2_FASTQ} --outdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/Quality/${LIBRARY_NAME}

printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### PhiX Alignment #######${NC}\n"
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}
bwa mem -k 16 -r 1 -M -T 15 -t ${MAXTHREADS} -v 1 ${PHIX_GENOME} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 25 -f 35 -uS - > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.PE.bam;
samtools view ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.PE.bam | cut -f 1 > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.list
sort --parallel=5 ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.list > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.sorted.list;
rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.list;
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
zcat ${R1_FASTQ} | python ${SCRIPT_DIR}/fqreverseextract.pureheader.py ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.sorted.list | pigz --best -f -c > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R1}_nophix.fastq.gz &
zcat ${R2_FASTQ} | python ${SCRIPT_DIR}/fqreverseextract.pureheader.py ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.sorted.list | pigz --best -f -c > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R2}_nophix.fastq.gz &
wait
NUMBER_PHIX_READS=`wc -l ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.sorted.list | cut -d' ' -f1 `;
printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${RED}##### PhiX READS: ${NUMBER_PHIX_READS} #####${NC}\n"
pigz -f ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/PhiX/${LIBRARY_NAME}/PhiX.header.sorted.list;
R1_FASTQ="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R1}_nophix.fastq.gz";
R2_FASTQ="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R2}_nophix.fastq.gz";


printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Ribosomal DNA Alignment #####${NC}\n"
mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}
bwa mem -k 16 -r 1 -M -T 15 -t ${MAXTHREADS} -v 1 ${RIBOSOMAL_GENOME_hum5SrDNA} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 25 -f 35 -uS - > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/hum5SrDNA.PE.bam &
bwa mem -k 16 -r 1 -M -T 15 -t ${MAXTHREADS} -v 1 ${RIBOSOMAL_GENOME_humRibosomal} <(zcat ${R1_FASTQ} ) <(zcat ${R2_FASTQ} ) | samtools view -F 2308 -q 25 -f 35 -uS - > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/humRibosomal.PE.bam &
wait
samtools view ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/hum5SrDNA.PE.bam | cut -f 1 > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/hum5SrDNA.list &
samtools view ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/humRibosomal.PE.bam | cut -f 1 > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/humRibosomal.list &
wait
cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/hum5SrDNA.list ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/humRibosomal.list > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.list;
sort --parallel=5 -u ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.list > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.sorted.list
rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/hum5SrDNA.list ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/humRibosomal.list ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.list;
BNAME_R1=`basename ${R1_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
BNAME_R2=`basename ${R2_FASTQ} | sed 's/.gz//g' | cut -d'.' -f1`;
zcat ${R1_FASTQ} | python ${SCRIPT_DIR}/fqreverseextract.pureheader.py ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.sorted.list | pigz --best -f -c > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R1}_noRibosomal.fastq.gz &
zcat ${R2_FASTQ} | python ${SCRIPT_DIR}/fqreverseextract.pureheader.py ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.sorted.list | pigz --best -f -c > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R2}_noRibosomal.fastq.gz &
wait
R1_FASTQ="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R1}_noRibosomal.fastq.gz";
R2_FASTQ="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/FastQ/${LIBRARY_NAME}/${BNAME_R2}_noRibosomal.fastq.gz";
NUMBER_RIBOSOMAL_READS=`wc -l ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.sorted.list | cut -d' ' -f1 `;
printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${RED}##### Ribosomal READS: ${NUMBER_RIBOSOMAL_READS} #####${NC}\n"
pigz -f ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/RibosomalRNA/${LIBRARY_NAME}/Ribosomal.header.sorted.list;


if [ ${ANALYSIS_PROTOCOL} = "tophat" ]; then
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC
	if [ ${LIBRARY_TYPE} = "fr-firststrand" ]; then
		printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Mapping Reads to the Genome #####${NC}\n"
		printf "${GREEN}@@@@ Splicing read mapping --> TopHat2${NC}\n"
		tophat2 -p ${MAXTHREADS} --library-type fr-firststrand -G ${GTF_FILE} --mate-inner-dist 0 --mate-std-dev 80 --no-coverage-search -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME} ${REFERENCE_GENOME_BOWTIE} ${R1_FASTQ} ${R2_FASTQ};
		samtools index ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/accepted_hits.bam;
		BAM="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/accepted_hits.bam";
	elif [ ${LIBRARY_TYPE} = "fr-secondstrand" ]; then
		printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Mapping Reads to the Genome #####${NC}\n"
		printf "${GREEN}@@@@ Splicing read mapping --> TopHat2${NC}\n"
		tophat2 -p ${MAXTHREADS} --library-type fr-secondstrand -G ${GTF_FILE} --mate-inner-dist 0 --mate-std-dev 80 --no-coverage-search -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME} ${REFERENCE_GENOME_BOWTIE} ${R1_FASTQ} ${R2_FASTQ};
		samtools index ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/accepted_hits.bam;
		BAM="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/accepted_hits.bam";
	fi
	READS_MAPPED=$((`cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/align_summary.txt | grep 'Mapped' | cut -d' ' -f17`)) ;
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${RED}##### READS MAPPED: ${READS_MAPPED} #####${NC}\n"
	cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv && echo "${LIBRARY_NAME},${BAM},${SAMPLE_TYPE}" >> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### RSeQC Report #######${NC}\n"
	inner_distance.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_annotation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_duplication.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_saturation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bam_stat.py -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.bam_stat.txt;
	read_distribution.py -r ${BED_FILE} -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.read_distribution.txt;
#	geneBody_coverage.py -r ${BED_FILE}  -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_quality.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bigwig -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/${LIBRARY_NAME}.bw
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bedgraph -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/TopHat2/${LIBRARY_NAME}/${LIBRARY_NAME}.bedgraph

elif [ ${ANALYSIS_PROTOCOL} = "hisat"  ]; then
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC
	if [ ${LIBRARY_TYPE} = "fr-firststrand" ]; then
		printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Mapping Reads to the Genome #####${NC}\n"
		printf "${GREEN}@@@@ Splicing read mapping --> HISAT2${NC}\n"
		hisat2 -p ${MAXTHREADS} --dta --new-summary --rna-strandness RF -x ${REFERENCE_GENOME_HISAT2} -1 ${R1_FASTQ} -2 ${R2_FASTQ} -S ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam 2> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.txt;
		samtools view -@ ${MAXTHREADS} -b -S ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
		rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam; 
		samtools sort -@ ${MAXTHREADS} ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
		rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
		samtools index ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
		BAM="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam";
	elif [ ${LIBRARY_TYPE} = "fr-secondstrand" ]; then
		printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Mapping Reads to the Genome #####${NC}\n"
		printf "${GREEN}@@@@ Splicing read mapping --> HISAT2${NC}\n"
		hisat2 -p ${MAXTHREADS} --dta --new-summary --rna-strandness FR -x ${REFERENCE_GENOME_HISAT2} -1 ${R1_FASTQ} -2 ${R2_FASTQ} -S ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam 2> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.txt;
		samtools view -@ ${MAXTHREADS} -b -S ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
		rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sam; 
		samtools sort -@ ${MAXTHREADS} ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
		rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
		samtools index ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
		BAM="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam";
	fi
	READS_MAPPED_1=$((`cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.txt | grep 'Aligned concordantly 1 time' | cut -d' ' -f5`)) ;
	READS_MAPPED_2=$((`cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.txt | grep 'Aligned concordantly >1 times' | cut -d' ' -f5`)) ;
	READS_MAPPED=$(($READS_MAPPED_1 + $READS_MAPPED_2));
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${RED}##### READS MAPPED: ${READS_MAPPED} #####${NC}\n"
	cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv && echo "${LIBRARY_NAME},${BAM},${SAMPLE_TYPE}" >> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### RSeQC Report #######${NC}\n"
	inner_distance.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_annotation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_duplication.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_saturation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bam_stat.py -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.bam_stat.txt;
	read_distribution.py -r ${BED_FILE} -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.read_distribution.txt;
#	geneBody_coverage.py -r ${BED_FILE}  -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_quality.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bigwig -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bw
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bedgraph -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/HISAT2/${LIBRARY_NAME}/${LIBRARY_NAME}.bedgraph

elif [ ${ANALYSIS_PROTOCOL} = "star"  ]; then
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}
	mkdir ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}##### Mapping Reads to the Genome #####${NC}\n"
	printf "${GREEN}@@@@ Splicing read mapping --> STAR${NC}\n"
	STAR --runMode alignReads  --runThreadN ${MAXTHREADS} --genomeDir ${REFERENCE_GENOME_STAR} --readFilesCommand zcat --outFileNamePrefix ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME} --readFilesIn ${R1_FASTQ} ${R2_FASTQ};
	samtools view -@ ${MAXTHREADS} -b -S ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}Aligned.out.sam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
	rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}Aligned.out.sam; 
	samtools sort -@ ${MAXTHREADS} ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.bam > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
	rm ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.bam;
	samtools index ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam;
	BAM="${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.sorted.bam";
	READS_MAPPED_1=$((`cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}Log.final.out | grep 'Uniquely mapped reads number' | cut -f2`)) ;
	READS_MAPPED_2=$((`cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}Log.final.out | grep 'Number of reads mapped to multiple loci' | cut -f2`)) ;
	READS_MAPPED=$(($READS_MAPPED_1 + $READS_MAPPED_2));
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${RED}##### READS MAPPED: ${READS_MAPPED} #####${NC}\n"
	cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv && echo "${LIBRARY_NAME},${BAM},${SAMPLE_TYPE}" >> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/input_all.csv
	printf "<`date +'%Y-%m-%d %H:%M:%S'`> ${YELLOW}####### RSeQC Report #######${NC}\n"
	inner_distance.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_annotation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_duplication.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	junction_saturation.py -r ${BED_FILE} -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bam_stat.py -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.bam_stat.txt;
	read_distribution.py -r ${BED_FILE} -i ${BAM} > ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME}.read_distribution.txt;
#	geneBody_coverage.py -r ${BED_FILE}  -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
#	read_quality.py -i ${BAM} -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/RSeQC/${LIBRARY_NAME};
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bigwig -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.bw
	bamCoverage -p ${MAXTHREADS} -b ${BAM} --outFileFormat bedgraph -o ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/STAR/${LIBRARY_NAME}/${LIBRARY_NAME}.bedgraph
fi

cat ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/reports/sample_report.csv && echo -e "${LIBRARY_NAME}\t${SAMPLE_TYPE}\t${NUMBER_RAW_READS}\t${NUMBER_PHIX_READS}\t${NUMBER_RIBOSOMAL_READS}\t${READS_MAPPED}" >> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}/reports/sample_report.csv

