
description=$(zenity --forms --title="RNA-Seq project" --text="Assign experiment code (no spaces!)" --add-entry="PROJECT NAME" --add-entry="POOL NAME" --add-entry="SAMPLE NAMES 
(first the controls,sep by ',')
[cntrl7,cntrl8,actd10,actd11]" --add-entry="SAMPLE TYPES (sep by ',')
[cntrl,cntrl,actd,actd]" --add-entry="LOG NAME")
pname=$(echo $description | cut -d'|' -f1)
poolname=$(echo $description | cut -d'|' -f2)
snames=$(echo $description | cut -d'|' -f3)
stype=$(echo $description | cut -d'|' -f4)
LOG=$(echo $description | cut -d'|' -f5)

array=$(echo $snames | tr "," "\n")
arr=($array)
READ1b=()
READ2b=()

se_pe=$(zenity --list --text="" --radiolist --column "" --column "" --hide-header --title="Paired end/Single end" TRUE "Paired_end" FALSE "Single_end")

for sample in `seq 1 "${#arr[@]}"`; do
	zenity --info --title="READ-1" --text="Select read-1 file for "${arr[$(($sample-1))]} --ok-label="OK";
	READ1a=$(zenity --file-selection --title="***READ-1***"  --text="Select read-1 file");
	READ1b+=( ${READ1a} )
done

if [ ${se_pe} = "Paired_end" ]; then
	for sample in `seq 1 "${#arr[@]}"`; do
		zenity --info --title="READ-2" --text="Select read-2 file for "${arr[$(($sample-1))]} --ok-label="OK";
		READ2a=$(zenity --file-selection --title="***READ-2***"  --text="Select read-2 file");
		READ2b+=( ${READ2a} )
	done
fi

function join { local IFS="$1"; shift; echo "$*"; }
READ1=$(join , ${READ1b[@]})
READ2=$(join , ${READ2b[@]})

zenity --info --title="BED file" --text="Select BED file" --ok-label="OK" 
BED=$(zenity --file-selection --filename /opt/genome/human/hg19/annotation/hg19.refseq.bed12 --title="***BED file***"  --text="Select BED file")

zenity --info --title="Phix genome" --text="Select Phix genome" --ok-label="OK" 
PHIX=$(zenity --file-selection --filename /opt/genome/control/phix174/bwa/phiX174.fa --title="***Phix genome***"  --text="Select Phix genome")

zenity --info --title="Ribosomal genome 1" --text="Select Ribosomal genome 1" --ok-label="OK" 
RIB1=$(zenity --file-selection --filename /opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa --title="***Ribosomal genome 1***"  --text="Select Ribosomal genome 1")

zenity --info --title="Ribosomal genome 2" --text="Select Ribosomal genome 2" --ok-label="OK" 
RIB2=$(zenity --file-selection --filename /opt/genome/human/hg19/contam/bwa/humRibosomal.fa --title="***Ribosomal genome 2***"  --text="Select Ribosomal genome 2")

zenity --info --title="GTF file" --text="Select GTF file" --ok-label="OK" 
GTF=$(zenity --file-selection --filename /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf --title="***GTF file***"  --text="Select GTF file")

zenity --info --title="Reference genome" --text="Select reference genome" --ok-label="OK" 
REF=$(zenity --file-selection --filename /opt/genome/human/hg19/index/hg19.fa --title="***Reference genome***"  --text="Select reference genome")

zenity --info --title="Script path" --text="Select script directory" --ok-label="OK" 
R=$(zenity --file-selection --directory --filename /opt/applications/src/arpir/ARPIR --title="***Script path***"  --text="Select script directory")

zenity --info --title="Output directory" --text="Select output directory" --ok-label="OK" 
OUT=$(zenity --file-selection --directory --title="***Output directory***"  --text="Select output directory")

zenity --info --title="Log directory" --text="Select log directory" --ok-label="OK" 
LOGF=$(zenity --file-selection --directory --filename /opt/ngs/logs --title="***Log directory***"  --text="Select log directory")

library=$(zenity --list --text="Choose library type" --radiolist --column "" --column "Library type" --hide-header --title="Library type" TRUE "fr-firststrand" FALSE "fr-secondstrand")

alignment=$(zenity --list --text="Choose alignment method" --radiolist --column "" --column "Alignment method" --hide-header --title="Alignment method" TRUE "hisat" FALSE "tophat")

zenity --info --title="Reference genome Bowtie" --text="Select reference genome (with indexes) for Bowtie" --ok-label="OK" 
REF_BOWTIE=$(zenity --file-selection --filename /opt/genome/human/hg19/index/bowtie2/hg19.fa --title="***Reference genome Bowtie***"  --text="Select reference genome for Bowtie")
REF_BOWTIE=$(echo $REF_BOWTIE | sed s/.fa//g)

zenity --info --title="Reference genome Hisat2" --text="Select reference genome (with indexes) for Hisat2" --ok-label="OK" 
REF_HISAT=$(zenity --file-selection --filename /opt/genome/human/hg19/index/hisat2/hg19.fa --title="***Reference genome Hisat2***"  --text="Select reference genome for Hisat2")
REF_HISAT=$(echo $REF_HISAT | sed s/.fa//g)

if [ ${alignment} = "hisat" ]; then
	quant="featureCounts";
	dea=$(zenity --list --text="Choose differential expression analysis method" --radiolist --column "" --column "DEA method" --hide-header --title="DEA method" TRUE "edgeR" FALSE "DESeq2");
fi

if [ ${alignment} = "tophat" ]; then
	quant="Cufflinks";
	dea="cummeRbund";
fi

meta=$(zenity --list --text="Choose whether to perform analysis with final meta-analysis or stop after quantification and DEA analysis" --radiolist --column "" --column "Analysis" --hide-header --title="Analysis" TRUE "full" FALSE "quant")

if [ ${meta} = "full" ]; then
	max_cat=$(zenity --forms --title="Max category" --text="Number of category to show in R plots" --add-entry="Category number");
elif [ ${meta} = "quant" ]; then
	max_cat="5"
fi

threads=$(zenity --forms --title="THREADS" --text="Number of threads" --add-entry="THREADS")

echo "project_name = "$pname"

pool_name = "$poolname"

sample_names = "$snames"

sample_types = "$stype"

log = "$LOG"

reads_1 = "$READ1"

reads_2 = "$READ2"

bowtie = "$REF_BOWTIE"

hisat = "$REF_HISAT"

BED_file = "$BED"

phix_genome = "$PHIX"

ribosomal_genomes = "$RIB1", "$RIB2"

GTF_file = "$GTF"

reference_genome = "$REF"

script_directory = "$R"

output_directory = "$OUT"

log_directory = "$LOGF"

meta-analysis = "$meta"

category_number = "$max_cat"

threads = "$threads | zenity --text-info --title="Summary" --width=700 --height=600 --ok-label="OK" --cancel-label="Cancel"

if [ ${se_pe} = "Paired_end" ]; then
	if [ "$?" -eq "0" ]; then
		cd ${OUT}
		python $R/ARPIR.py -n $pname -pn $poolname -sn $snames -r1 $READ1 -r2 $READ2 -type $stype -rb $REF_BOWTIE -rh $REF_HISAT -bed $BED -ph $PHIX -rib1 $RIB1 -rib2 $RIB2 -t $threads -g $GTF -a $alignment -l $library -q $quant -r $REF -dea $dea -script_path $R -o $OUT -meta $meta -cat $max_cat 2>&1 >> ${LOGF}/${LOG}.log
	fi
fi

if [ ${se_pe} = "Single_end" ]; then
	if [ "$?" -eq "0" ]; then
		cd ${OUT}
		python $R/ARPIR.py -n $pname -pn $poolname -sn $snames -r1 $READ1 -type $stype -rb $REF_BOWTIE -rh $REF_HISAT -bed $BED -ph $PHIX -rib1 $RIB1 -rib2 $RIB2 -t $threads -g $GTF -a $alignment -l $library -q $quant -r $REF -dea $dea -script_path $R -o $OUT -meta $meta -cat $max_cat 2>&1 >> ${LOGF}/${LOG}.log
	fi
fi

