#!/bin/bash

cd ../
mkdir TopHat
mkdir Hisat
mkdir Hisat/StringTie
mkdir kallisto

## TOPHAT ##
touch TopHat/cufflinks.transcripts.allsamples.txt
for file in *_1.fastq.gz; do
tophat2 -p 10 --library-type fr-firststrand -G /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf --mate-inner-dist 0 --mate-std-dev 80 --no-coverage-search -o TopHat/"$(basename "$file" _1.fastq.gz)" /opt/genome/human/hg19/index/bowtie2/hg19 "$file" "$(basename "$file" 1.fastq.gz)2.fastq.gz";
samtools index TopHat/"$(basename "$file" _1.fastq.gz)"/accepted_hits.bam; 

## CUFFLINKS ## 
cufflinks --label "$(basename "$file" _1.fastq.gz)" --no-update-check --library-type fr-firststrand --verbose --GTF /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf -p 10 -o TopHat/"$(basename "$file" _1.fastq.gz)"/cufflinks TopHat/"$(basename "$file" _1.fastq.gz)"/accepted_hits.bam; 
cat TopHat/cufflinks.transcripts.allsamples.txt && echo TopHat/"$(basename "$file" _1.fastq.gz)"/cufflinks/transcripts.gtf >> TopHat/cufflinks.transcripts.allsamples.txt; done
cuffmerge -g /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf -s /opt/genome/human/hg19/index/hg19.fa -p 10 -o TopHat TopHat/cufflinks.transcripts.allsamples.txt
for file in *_1.fastq.gz; do
cuffquant -p10 --library-type  fr-firststrand --no-update-check -o TopHat/"$(basename "$file" _1.fastq.gz)"/cuffquant TopHat/merged.gtf TopHat/"$(basename "$file" _1.fastq.gz)"/accepted_hits.bam; done
cuffdiff -p 10 -L cntrl,treat -o TopHat/cuffdiff --library-type  fr-firststrand --no-update-check TopHat/merged.gtf TopHat/sample_01/cuffquant/abundances.cxb,TopHat/sample_02/cuffquant/abundances.cxb,TopHat/sample_03/cuffquant/abundances.cxb TopHat/sample_04/cuffquant/abundances.cxb,TopHat/sample_05/cuffquant/abundances.cxb,TopHat/sample_06/cuffquant/abundances.cxb


## HISAT ##
for file in *_1.fastq.gz; do
hisat2 -p 10 --dta --rna-strandness RF -x /opt/genome/human/hg19/index/hisat2/hg19 -1 "$file" -2 "$(basename "$file" 1.fastq.gz)2.fastq.gz" -S Hisat/"$(basename "$file" _1.fastq.gz).sam";
samtools view -@ 10 -b -S Hisat/"$(basename "$file" _1.fastq.gz).sam" > Hisat/"$(basename "$file" _1.fastq.gz).bam";
samtools sort -@ 10 Hisat/"$(basename "$file" _1.fastq.gz).bam" > Hisat/"$(basename "$file" _1.fastq.gz).sorted.bam";
samtools index Hisat/"$(basename "$file" _1.fastq.gz).sorted.bam"; done

## FEATURECOUNTS ##
featureCounts -T 10 -p -a /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf -o Hisat/featureCounts_counts.txt Hisat/sample_01.sorted.bam Hisat/sample_02.sorted.bam Hisat/sample_03.sorted.bam Hisat/sample_04.sorted.bam Hisat/sample_05.sorted.bam Hisat/sample_06.sorted.bam 

## STRINGTIE ##
for file in *_1.fastq.gz; do
stringtie -p 10 -G /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf -e -B -o Hisat/StringTie/"$(basename "$file" _1.fastq.gz)"/transcripts.gtf -A Hisat/StringTie/"$(basename "$file" _1.fastq.gz)"/gene_abundances.tsv Hisat/"$(basename "$file" _1.fastq.gz).sorted.bam";
perl Scripts/stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs=Hisat/StringTie/"$(basename "$file" _1.fastq.gz)" --transcript_matrix_file=Hisat/StringTie/"$(basename "$file" _1.fastq.gz)"/transcript_fpkm.tsv --gene_matrix_file=Hisat/StringTie/"$(basename "$file" _1.fastq.gz)"/gene_fpkm.tsv; done


## KALLISTO ##
for file in *_1.fastq.gz; do
kallisto quant -i /opt/genome/human/hg19/index/kallisto/hsGRCh37.75_kallisto -t 10 -b 100 "$file" "$(basename "$file" 1.fastq.gz)2.fastq.gz" -o kallisto/"$(basename "$file" _1.fastq.gz)"; done
cd kallisto
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv
