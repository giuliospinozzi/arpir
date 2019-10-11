#!/bin/bash

## Space
	LEFTSPACE=`df -P | grep 'sda' | awk '{print $4}'`
	if [ ${LEFTSPACE} -lt 25600000 ]; then 
		echo "
    
        *****************************************************************
        |                                                               |
        |   YOU DO NOT HAVE ENOUGH SPACE LEFT ON HARD DISK!!            |
        |                                                               |
        |   I WILL NOT PROCEED...                                       |
        |                                                               |       
        |   FREE SPACE BEFORE, LEAVING AT LEAST 25GB ON HD              |
        |                                                               |
        *****************************************************************
    
        ";
        exit
    fi


######################  COLOR  ######################

	RED='\033[0;31m' 
	NC='\033[0m' # No Color
	GREEN='\033[0;32m'
	YELLOW='\033[1;33m'
	CYAN='\033[0;36m'
	BLUE='\033[0;34m'
	PURPLE='\033[0;35m'

######################  MY FUNCTION  ######################

	function installDepend () {
		apt-get install -y $1
	}

	function installPyPack () {
		pip install $1
	}

	function checkAllPacks () {

		echo "@  Checking if all pack are installed"
		list_dependencies=(zenity libcurl4-openssl-dev libmagick++-dev libmariadbclient-dev libssl-dev pandoc fastqc python-pip bwa samtools tophat hisat2 cufflinks)
		for i in ${list_dependencies[@]}; do
			if [[ `dpkg -s $i 2> /dev/null | wc -l` > 0 ]]; then
				printf "[ ${GREEN}OK${NC} ] $i installed correctly\n"
			else
				printf "[ ${RED}MISSING${NC} ] $i not installed correctly\n"
				installDepend "$i"
			fi
		done
	}

######################  HELP  ######################

	usage()
	{
		echo
		echo "Setting up ARPIR"
		echo
		echo "#Summary:"
		echo "	This program is an introduction for ARPIR"
		echo "	You will be able to set up your machine to allow work ARPIR correctly"
		echo "	Will downloaded all the dependencies required; $0 installs packages to their own directory"
		echo "	If you prefer to configure manually your machine visit the link: https://github.com/giuliospinozzi/arpir"
		echo
		echo "#Usage: sudo bash conf.sh"
		echo
		exit 
	}

######################  MAIN  ######################

###################### DEPENDENCIES ######################

	list_dependencies=(zenity r-base libcurl4-openssl-dev libmagick++-dev libmariadbclient-dev libssl-dev pandoc fastqc python-pip bwa samtools tophat hisat2 cufflinks pigz libgd-dev)
	list_dependencies1=(argparse pandas multiqc RSeQC)

	echo "STEP 1: DEPENDENCIES"
	echo "Installing missing dependencies"
	echo
	
	echo "Updating your system..."
	apt-get update -y
	apt-get upgrade -y

	for i in "${list_dependencies[@]}"; do 
		
		if [[ $i == "r-base" ]]; then

			if [[( `dpkg -s r-base 2> /dev/null | wc -l` > 0 || -e /usr/bin/R )]]; then 
				printf "[ ${GREEN}OK${NC} ] R already exists\n"
				continue
			else
				printf "[ ${RED}NO${NC} ] R not found, Installing ${i}\n"
				apt install r-base -y 2> /dev/null
			fi	
		else
			if [[ `dpkg -s $i 2> /dev/null | wc -l ` > 0 ]]; then 
				printf "[ ${GREEN}OK${NC} ] $i already exists\n"
			else
				printf "[ ${RED}NO${NC} ] $i not found, Installing ${i}\n"
				installDepend "$i"
			fi
		fi
	done

	echo
	checkAllPacks

	for i in "${list_dependencies1[@]}"; do
		printf " Installing ${i}\n" 
		installPyPack "$i"
	done

	curl -OL "https://github.com/giuliospinozzi/arpir/raw/master/genome.zip"
	unzip genome.zip
	rm genome.zip

	if [[( -e /usr/bin/fastq_screen || -e /usr/local/bin/fastq_screen )]]; then
		printf "[ ${GREEN}OK${NC} ] FastQ Screen already exists\n"
		continue
	else
		printf "[ ${RED}NO${NC} ] FastQ Screen not found, Installing FastQ Screen\n"
		curl -OL "https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/fastq_screen_v0.14.0.tar.gz"
		tar -xf fastq_screen_v0.14.0.tar.gz
		rm fastq_screen_v0.14.0.tar.gz
		mv fastq_screen_v0.14.0 fastq_screen
		mv fastq_screen /etc 2> /dev/null
		chmod +x /etc/fastq_screen/fastq_screen
		ln -sf /etc/fastq_screen/fastq_screen /usr/bin/fastq_screen 2> /dev/null
		mv genome/fastq_screen.conf /etc/fastq_screen
	fi

	if [[( -e /usr/bin/featureCounts || -e /usr/local/bin/featureCounts )]]; then
		printf "[ ${GREEN}OK${NC} ] featureCounts already exists\n"
		continue
	else
		printf "[ ${RED}NO${NC} ] featureCounts not found, Installing featureCounts\n"
		curl -OL "https://sourceforge.net/projects/subread/files/subread-1.5.3/subread-1.5.3-Linux-x86_64.tar.gz"
		tar -xf subread-1.5.3-Linux-x86_64.tar.gz
		rm subread-1.5.3-Linux-x86_64.tar.gz
		mv subread-1.5.3-Linux-x86_64 subread
		mv subread /etc 2> /dev/null
		chmod +x /etc/subread/bin/featureCounts
		ln -sf /etc/subread/bin/featureCounts /usr/bin/featureCounts 2> /dev/null
	fi

	printf "Installing R packages\n"
	R -e 'install.packages(c("ggfortify", "ggrepel", "RColorBrewer", "gplots", "dplyr", "igraph", "scales", "treemap", "shiny", "DT", "magick", "rlist", "visNetwork", "shinyjs", "knitr"));install.packages("https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.0.9.tar.gz", repos=NULL, type="source");source("https://bioconductor.org/biocLite.R");biocLite(c("cummeRbund", "edgeR", "DESeq2", "genefilter", "clusterProfiler", "org.Hs.eg.db", "pathview"))'
	yes | perl -MCPAN -e "install GD"
	yes | perl -MCPAN -e "install GD::Graph"


###################### Download Genome ######################

	echo ""
	echo "STEP 2: Download Genome"
	echo ""

	#bed and gtf files
	echo "Downloading annotation files..."
	mkdir -p /opt/genome/human/hg19/annotation
	mv genome/hg19* /opt/genome/human/hg19/annotation

	# phix genome
	echo "Downloading PhiX genome..."
	mkdir -p /opt/genome/control/phix174/bwa/
	mv genome/phiX174.fa /opt/genome/control/phix174/bwa
	bwa index /opt/genome/control/phix174/bwa/phiX174.fa

	# ribosomal genome
	echo "Downloading ribosomal genomes..."
	mkdir -p /opt/genome/human/hg19/contam/bwa/
	mv genome/hum* /opt/genome/human/hg19/contam/bwa
	bwa index /opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa
	bwa index /opt/genome/human/hg19/contam/bwa/humRibosomal.fa
	rm -r genome

	# reference genome
	echo "Downloading reference genome (hg19)..."
	mkdir -p /opt/genome/human/hg19/index/bwa/
	cd /opt/genome/human/hg19/index
	curl -OL "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
	gunzip hg19.fa.gz
	samtools faidx /opt/genome/human/hg19/index/hg19.fa
	ln -s /opt/genome/human/hg19/index/hg19.fa /opt/genome/human/hg19/index/bwa/hg19.fa
	bwa index /opt/genome/human/hg19/index/bwa/hg19.fa

	# bowtie index
	echo "Creating bowtie2 index..."
	mkdir -p /opt/genome/human/hg19/index/bowtie2
	ln -s /opt/genome/human/hg19/index/hg19.fa /opt/genome/human/hg19/index/bowtie2/hg19.fa
	bowtie2-build /opt/genome/human/hg19/index/bowtie2/hg19.fa /opt/genome/human/hg19/index/bowtie2/hg19

	#hisat index
	echo "Creating hisat2 index..."
	mkdir -p /opt/genome/human/hg19/index/hisat2
	ln -s /opt/genome/human/hg19/index/hg19.fa /opt/genome/human/hg19/index/hisat2/hg19.fa
	hisat2-build /opt/genome/human/hg19/index/hisat2/hg19.fa /opt/genome/human/hg19/index/hisat2/hg19

	# # mouse genome
	# sudo mkdir -p /opt/genome/mouse/mm10/index/bwa/
	# sudo cd /opt/genome/mouse/mm10/index/bwa/
	# sudo curl -OL "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz"
	# sudo tar -xvzf chromFa.tar.gz
	# sudo rm chromFa.tar.gz 
	# for i in $(seq 1 19); do
	# 	cat chr$i.fa >> mm10.fa
	# done
	# cat chrX.fa >> mm10.fa
	# cat chrY.fa >> mm10.fa
	# cat chrM.fa >> mm10.fa
	# sudo rm chr*
	# bwa index /opt/genome/mouse/mm10/index/bwa/mm10.fa







