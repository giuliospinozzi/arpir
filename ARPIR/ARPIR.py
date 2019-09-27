#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, argparse, sys, csv
import pandas as pd

header = """
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
  |  Version:  1.0 - CREO - RNAseq                         |
  +--------------------------------------------------------+

  optional arguments (NO SPACES!!):
  
  -n, --project-name
  -pn, --pool-name
  -sn, --sample-name
  -r1, --read1
  -r2, --read2
  -type, --type [cntrl,treat]
  -rb, --reference-genome-bowtie [/opt/genome/human/hg19/index/bowtie2/hg19]
  -rh, --reference-genome-hisat2 [/opt/genome/human/hg19/index/hisat2/hg19]
  -bed, --bed-file [/opt/genome/human/hg19/annotation/hg19.refseq.bed12]
  -ph, --phix-genome [/opt/genome/control/phix174/bwa/phiX174.fa]
  -rib1, --ribosomal-genome1 [/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa]
  -rib2, --ribosomal-genome2 [/opt/genome/human/hg19/contam/bwa/humRibosomal.fa]
  -t, --threads [12]
  -g, --gtf [/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf]
  -a, --alignment-method [hisat,tophat]
  -l, --library-type [fr-firststrand,fr-secondstrand]
  -q, --quantification-method [featureCounts,Cufflinks]
  -r, --reference-genome [/opt/genome/human/hg19/index/hg19.fa]
  -dea, --dea-method [edgeR,DESeq2,cummeRbund]
  -r_path, --r_path [/opt/applications/src/arpir/ARPIR]
  -o, --output-dir
  -meta, --meta-analysis [full,quant]
  -cat, --max-cat [5]
  -comp, --comparisons
  
""" 


description = "This application makes quality control, pre-processing, alignment, transcript quantification, differential expression analysis and optionaly meta-analysis on FastQ files"

usage_example = """
Examples of usage:
 (1) RNA-Seq analysis: 
    <appname> -n <project_name> -pn <pool_name> -sn <sample_A,sample_B,sample_C,sample_D> -r1 <sample_A_read1,sample_B_read1,sample_C_read1,sample_D_read1> -r2 <sample_A_read2,sample_B_read2,sample_C_read2,sample_D_read2> -type <cntrl,cntrl,treat,treat> -o <output_directory> [options]*

"""

print header, description, usage_example

parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ RNA-Seq analysis: from BAM to DEA ] \n", description = description)
parser.add_argument('-n', '--project-name', dest="project_name", help="Project name. \n No default option. \n", action="store", required=True)
parser.add_argument('-pn', '--pool-name', dest="pool_name", help="Pool name. \n No default option. \n", action="store", required=True)
parser.add_argument('-sn', '--sample-name', dest="sample_name", help="Sample names (',' sep, first the control samples). \n No default option. \n", action="store", required=True)
parser.add_argument('-r1', '--read1', dest="read1", help="Read 1 fastq path (',' sep). \n No default option. \n", action="store", required=True)
parser.add_argument('-r2', '--read2', dest="read2", help="Read 2 fastq path (',' sep). \n Default: single-end. \n", action="store", required=False, default="")
parser.add_argument('-type', '--type', dest="stype", help="Sample types (cntrl,treat). \n No default option. \n", action="store", required=True)
parser.add_argument('-rb', '--reference-genome-bowtie', dest="ref_bowtie", help="Reference genome file path for bowtie. \n Default: /opt/genome/human/hg19/index/bowtie2/hg19. \n", action="store", required=False, default="/opt/genome/human/hg19/index/bowtie2/hg19")
parser.add_argument('-rh', '--reference-genome-hisat2', dest="ref_hisat2", help="Reference genome file path for hisat2. \n Default: /opt/genome/human/hg19/index/hisat2/hg19. \n", action="store", required=False, default="/opt/genome/human/hg19/index/hisat2/hg19")
parser.add_argument('-bed', '--bed-file', dest="bed_file", help="Reference genome annotation file path. \n Default: /opt/genome/human/hg19/annotation/hg19.refseq.bed12. \n", action="store", required=False, default="/opt/genome/human/hg19/annotation/hg19.refseq.bed12")
parser.add_argument('-ph', '--phix-genome', dest="phix", help="Phix genome file path. \n Default: /opt/genome/control/phix174/bwa/phiX174.fa. \n", action="store", required=False, default="/opt/genome/control/phix174/bwa/phiX174.fa")
parser.add_argument('-rib1', '--ribosomal-genome1', dest="rib1", help="Ribosomal genome file path. \n Default: /opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa. \n", action="store", required=False, default="/opt/genome/human/hg19/contam/bwa/hum5SrDNA.fa")
parser.add_argument('-rib2', '--ribosomal-genome2', dest="rib2", help="Ribosomal genome file path. \n Default: /opt/genome/human/hg19/contam/bwa/humRibosomal.fa. \n", action="store", required=False, default="/opt/genome/human/hg19/contam/bwa/humRibosomal.fa")
parser.add_argument('-t', '--threads', dest="Threads", help="Max thread number. \n Default: 12. \n", action="store", required=False, default=12)
parser.add_argument('-g', '--gtf', dest="GTF", help="GTF file path. \n Default: /opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf. \n", action="store", required=False, default="/opt/genome/human/hg19/annotation/hg19.refgene.sorted.gtf")
parser.add_argument('-a', '--alignment-method', dest="a_method", help="Alignment method. \n Default: hisat; alternative: tophat. \n", action="store", required=False, default="hisat")
parser.add_argument('-l', '--library-type', dest="library_type", help="Library type. \n Default: fr-firststrand; alternative: fr-secondstrand. \n", action="store", required=False, default="fr-firststrand")
#parser.add_argument('-i', '--input', dest="input_path", help="Input path dataframe. \n No default option. \n", action="store", required=True)
parser.add_argument('-q', '--quantification-method', dest="q_method", help="Quantification method. \n Default: featureCounts; alternative: Cufflinks. \n", action="store", required=False, default="featureCounts")
parser.add_argument('-r', '--reference-genome', dest="ref_gen", help="Reference genome file path (only for Cufflinks). \n Default: human hg19. \n", action="store", required=False, default="/opt/genome/human/hg19/index/hg19.fa")
parser.add_argument('-dea', '--dea-method', dest="dea_method", help="Differential Expression Analysis method. \n Default: edgeR; alternatives: DESeq2, cummeRbund. \n", action="store", required=False, default="edgeR")
parser.add_argument('-r_path', '--r_path', dest="R_path", help="Script directory (alignment, quantification and DEA). \n Default: /opt/applications/src/arpir/ARPIR. \n", action="store", required=False, default="/opt/applications/src/arpir/ARPIR")
parser.add_argument('-o', '--output-dir', dest="output_dir", help="Output directory. \n No default option. \n", action="store", required=True)
parser.add_argument('-meta', '--meta-analysis', dest="meta", help="Analysis with or without final meta-analysis. \n Default: full; alternative: quant. \n", action="store", required=False, default="full")
parser.add_argument('-cat', '--max_cat', dest="max_cat", help="Max number of category showed in R plots for meta-analysis. \n Default: 5. \n", action="store", required=False, default="5")
parser.add_argument('-comp', '--comparisons', dest="comp", help="Comparisons (cntrl_VS_treat1,cntrl_VS_treat2). \n No default option. \n", action="store", required=True)

args = parser.parse_args()


########################################################################
####### GLOBAL VARS
########################################################################

########################################################################
####### FUNCTIONS
########################################################################

def checkArgs(args):
    """
    Check file path
    """
    print "Project name =", args.project_name, "\n"
    print "Pool names =", args.pool_name, "\n"
    print "Sample name =", args.sample_name, "\n"
    print "Read 1 fastq =", args.read1, "\n"
    print "Read 2 fastq =", args.read2, "\n"
    print "Sample types =", args.stype, "\n"
    print "Reference genome for bowtie =", args.ref_bowtie, "\n"
    print "Reference genome for hisat2 =", args.ref_hisat2, "\n"
    print "Bed file =", args.bed_file, "\n"
    print "Phix genome file =", args.phix, "\n"
    print "Ribosomal genome file 1 =", args.rib1, "\n"
    print "Ribosomal genome file 2 =", args.rib2, "\n"
    print "Max threads =", args.Threads, "\n"
    print "GTF file =", args.GTF, "\n"
    print "Alignment method =", args.a_method, "\n"
    print "Library type =", args.library_type, "\n"
    print "Quantification method =", args.q_method, "\n"
    print "Reference genome =", args.ref_gen, "\n"
    print "DEA method =", args.dea_method, "\n"
    print "Script path =", args.R_path, "\n"
    print "Output directory =", args.output_dir, "\n"
    print "Comparisons =", args.comp, "\n"
    if not os.path.isfile(args.bed_file):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid path for BED file.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(args.phix):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid path for Phix genome.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(args.rib1) or not os.path.isfile(args.rib2):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid paths for ribosomal genome.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(args.GTF):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid path for GTF.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isdir(args.output_dir):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid path for output dir.\n\tExit\n" + '\033[0m')
        sys.exit()
    if args.a_method == "hisat":
        h=args.ref_hisat2+".fa"
        if not os.path.isfile(h):
            print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid paths for index genome (hisat).\n\tExit\n" + '\033[0m')
            sys.exit()
    if args.a_method == "tophat":
        b=args.ref_bowtie+".fa"
        if not os.path.isfile(b):
            print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid paths for index genome (bowtie).\n\tExit\n" + '\033[0m')
            sys.exit()
    r1=args.read1.split(",")
    for i in range(0, (len(r1))):
        if not os.path.isfile(r1[i]):
            print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid paths for read 1.\n\tExit\n" + '\033[0m')
            sys.exit()
    if args.read2 != "":
        r2=args.read2.split(",")
        for j in range(0, (len(r2))):
            if not os.path.isfile(r2[j]):
                print ('\033[0;31m' + "\n[AP]\tError while reading files: no valid paths for read 2.\n\tExit\n" + '\033[0m')
                sys.exit()
    script_path1=os.path.abspath(os.path.dirname(sys.argv[0]))+"/alignment.sh"
    script_path2=os.path.abspath(os.path.dirname(sys.argv[0]))+"/quantification.py"
    script_path3=os.path.abspath(os.path.dirname(sys.argv[0]))+"/GO_pathway.R"
    script_path4=os.path.abspath(os.path.dirname(sys.argv[0]))+"/CNETPLOT_FUNCTION.R"
    script_path5=os.path.abspath(os.path.dirname(sys.argv[0]))+"/alignment_se.sh"
    script_path6=os.path.abspath(os.path.dirname(sys.argv[0]))+"/runcummeRbund.R"
    script_path7=os.path.abspath(os.path.dirname(sys.argv[0]))+"/runDESeq2.R"
    script_path8=os.path.abspath(os.path.dirname(sys.argv[0]))+"/runedgeR.R"
    script_path9=os.path.abspath(os.path.dirname(sys.argv[0]))+"/fqreverseextract.pureheader.py"
    if not os.path.isfile(script_path1) or not os.path.isfile(script_path5) or not os.path.isfile(script_path9):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no alignment scripts.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(script_path2):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no quantification script.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(script_path3) or not os.path.isfile(script_path4):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no meta-analysis scripts.\n\tExit\n" + '\033[0m')
        sys.exit()
    if not os.path.isfile(script_path6) or not os.path.isfile(script_path7) or not os.path.isfile(script_path8):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: no DEA scripts.\n\tExit\n" + '\033[0m')
        sys.exit()


def checkOptions(a_method,q_method,dea_method):
    """
    Check quantification and DEA method options
    """
    if a_method == 'hisat' and q_method == 'Cufflinks':
        print ('\033[0;31m' + "\n[AP]\tError: alignment method HISAT2 requires quantification method featureCounts\n\tExit\n" + '\033[0m')
        sys.exit()
    if a_method == 'tophat' and q_method == 'featureCounts':
        print ('\033[0;31m' + "\n[AP]\tError: alignment method Tophat requires quantification method Cufflinks\n\tExit\n" + '\033[0m')
        sys.exit()
    if q_method == 'Cufflinks' and dea_method == 'edgeR':
        print ('\033[0;31m' + "\n[AP]\tError: quantification method Cufflinks requires DEA method cummeRbund\n\tExit\n" + '\033[0m')
        sys.exit()
    if q_method == 'Cufflinks' and dea_method == 'DESeq2':
        print ('\033[0;31m' + "\n[AP]\tError: quantification method Cufflinks requires DEA method cummeRbund\n\tExit\n" + '\033[0m')
        sys.exit()
    if q_method == 'featureCounts' and dea_method == 'cummeRbund':
        print ('\033[0;31m' + "\n[AP]\tError: quantification method featureCounts requires DEA method edgeR or DESeq2\n\tExit\n" + '\033[0m')
        sys.exit()


def checkFile(read1,read2,stype,sample_name):
    """
    Check read and sample type number and case/control
    """
    read1a=read1.split(",")
    stype1=stype.split(",")
    sample_name1=sample_name.split(",")
    if len(read1a)!=len(stype1) or len(read1a)!=len(sample_name1):
        print ('\033[0;31m' + "\n[AP]\tError while reading files: read1, type and sample name must be of the same length.\n\tExit\n" + '\033[0m')
        sys.exit()
    if read2 != "":
        read2a=read2.split(",")
        if len(read1a)!=len(read2a):
            print ('\033[0;31m' + "\n[AP]\tError while reading files: read1 and read2 must be of the same length.\n\tExit\n" + '\033[0m')
            sys.exit()
    if len(set(stype1))<2:
        print ('\033[0;31m' + "\n[AP]\tError while reading files: there must be at least one case and one control.\n\tExit\n" + '\033[0m')
        sys.exit()
    


def alignment(project_name,pool_name,sample_name,output_dir,read1,read2,Threads,ref_bowtie,ref_hisat2,bed_file,phix,rib1,rib2,a_method,GTF,library_type,R_path,stype,q_method,dea_method,comp):
    """
    Run alignment script as desired
    """
#    print ('\033[1;33m' + "\n^^^^^^^^^Alignment script running^^^^^^^^^\n" + '\033[0m')
    cmd1="mkdir "+output_dir+"/"+project_name
    os.system(cmd1)
    cmd2="mkdir "+output_dir+"/"+project_name+"/"+pool_name
    os.system(cmd2)
    cmd3="mkdir "+output_dir+"/"+project_name+"/"+pool_name+"/reports"
    os.system(cmd3)    
    with open(output_dir+'/'+project_name+'/'+pool_name+'/input_all.csv', 'wb') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',')
        filewriter.writerow(['sample_name','BAM_path','Type'])
        csvfile.close()
    with open(output_dir+'/'+project_name+'/'+pool_name+'/reports/general_report.csv', 'wb') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        filewriter.writerow(['project_name',project_name])
        filewriter.writerow(['pool_name',pool_name])
        filewriter.writerow(['alignment_method',a_method])
        filewriter.writerow(['library_type',library_type])
        filewriter.writerow(['quantification_method',q_method])
        filewriter.writerow(['dea_method',dea_method])
        filewriter.writerow(['comparisons',comp])
        csvfile.close()
    if read2 == "":
        with open(output_dir+'/'+project_name+'/'+pool_name+'/reports/general_report.csv', 'a') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t')
            filewriter.writerow(['Single-end/Paired-end',"Single_end"])
            csvfile.close()
    if read2 != "":
        with open(output_dir+'/'+project_name+'/'+pool_name+'/reports/general_report.csv', 'a') as csvfile:
            filewriter = csv.writer(csvfile, delimiter='\t')
            filewriter.writerow(['Single-end/Paired-end',"Paired_end"])
            csvfile.close()
    with open(output_dir+'/'+project_name+'/'+pool_name+'/reports/sample_report.csv', 'wb') as csvfile:
        filewriter = csv.writer(csvfile, delimiter='\t')
        filewriter.writerow(['sample_name','sample_type','number_raw_reads','number_phix_reads','number_ribosomal_reads','number_mapped_reads'])
        csvfile.close()
    read1a=read1.split(",")
    sample_name1=sample_name.split(",")
    stype1=stype.split(",")
    if read2 != "":
        read2a=read2.split(",")
        for i in range(0, (len(read1a))):
            cmd="bash "+R_path+"/alignment.sh "+project_name+" "+pool_name+" "+sample_name1[i]+" "+output_dir+" "+read1a[i]+" "+read2a[i]+" "+str(Threads)+" "+ref_bowtie+" "+ref_hisat2+" "+bed_file+" "+phix+" "+rib1+" "+rib2+" "+a_method+" "+GTF+" "+library_type+" "+stype1[i]+" "+R_path
            os.system(cmd)
    if read2 == "":
        for i in range(0, (len(read1a))):
            cmd="bash "+R_path+"/alignment_se.sh "+project_name+" "+pool_name+" "+sample_name1[i]+" "+output_dir+" "+read1a[i]+" "+str(Threads)+" "+ref_bowtie+" "+ref_hisat2+" "+bed_file+" "+phix+" "+rib1+" "+rib2+" "+a_method+" "+GTF+" "+library_type+" "+stype1[i]+" "+R_path
            os.system(cmd)


def quantification(output_dir,project_name,pool_name,R_path,dea_method,q_method,Threads,GTF,library_type,ref_gen,comp):
    """
    Run quantification script as desired
    """
#    print ('\033[1;33m' + "\n^^^^^^^^^Quantification script running^^^^^^^^^\n" + '\033[0m')
    res_dir=output_dir+"/"+project_name+"/"+pool_name+"/"+comp
    cmd="python "+R_path+"/quantification.py -i "+output_dir+'/'+project_name+'/'+pool_name+'/input.csv'+" -o "+res_dir+"/Quantification_and_DEA"+" -r_path "+R_path+" -dea "+dea_method+" -q "+q_method+" -t "+str(Threads)+" -g "+GTF+" -l "+library_type+" -r "+ref_gen+" -comp "+comp
    os.system(cmd)
    
    
def rminput(output_dir,project_name,pool_name):
    """
    Remove input file
    """
    os.system("rm "+output_dir+'/'+project_name+'/'+pool_name+'/input.csv')
    
    
def metaanalysis(output_dir,R_path,project_name,pool_name,dea_method,max_cat,comp):
    """
    Run meta-analysis script as desired
    """
    cmd="Rscript --vanilla --verbose "+R_path+"/GO_pathway.R "+output_dir+'/'+project_name+'/'+pool_name+'/'+comp+'/Quantification_and_DEA/*-results.csv '+dea_method+" "+output_dir+'/'+project_name+'/'+pool_name+'/'+comp+"/Meta-analysis "+max_cat+" "+R_path
    os.system(cmd)


def report(output_dir,project_name,pool_name):
    """
    Generate multiqc report
    """
    cmd="multiqc --exclude bowtie2 -o "+output_dir+'/'+project_name+'/'+pool_name+"/reports/ "+output_dir+'/'+project_name+'/'+pool_name+"/"
    os.system(cmd)
        

#########################################################################################
### MAIN
#########################################################################################

def main():
    """
    Main part of the program.
    """
    # first check args and file paths
    checkArgs(args)
    checkOptions(args.a_method,args.q_method,args.dea_method)
    checkFile(args.read1,args.read2,args.stype,args.sample_name)

    # alignment
    alignment(args.project_name,args.pool_name,args.sample_name,args.output_dir,args.read1,args.read2,args.Threads,args.ref_bowtie,args.ref_hisat2,args.bed_file,args.phix,args.rib1,args.rib2,args.a_method,args.GTF,args.library_type,args.R_path,args.stype,args.q_method,args.dea_method,args.comp)

    comp1=args.comp.split(",")
    for i in range(0, (len(comp1))):
        cat_type=comp1[i].split("_VS_")
        data=pd.DataFrame.from_csv(args.output_dir+'/'+args.project_name+'/'+args.pool_name+'/input_all.csv',sep=',',index_col=None)
        input_all=data.loc[(data['Type'] == cat_type[0]) | (data['Type'] == cat_type[1])]
        input_all.to_csv (r''+args.output_dir+'/'+args.project_name+'/'+args.pool_name+'/input.csv', index = None, header=True)
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i])
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i]+"/Quantification_and_DEA")
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i]+"/Meta-analysis")
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i]+"/Meta-analysis/Gene_ontology")
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i]+"/Meta-analysis/Pathway_analysis")
        os.system("mkdir "+args.output_dir+"/"+args.project_name+"/"+args.pool_name+"/"+comp1[i]+"/Meta-analysis/Pathway_analysis/pathview")

        # quantification
        quantification(args.output_dir,args.project_name,args.pool_name,args.R_path,args.dea_method,args.q_method,args.Threads,args.GTF,args.library_type,args.ref_gen,comp1[i])

        # remove input file
        rminput(args.output_dir,args.project_name,args.pool_name)

        # meta-analysis (GO and pathway analysis)
        if args.meta == 'full':
            metaanalysis(args.output_dir,args.R_path,args.project_name,args.pool_name,args.dea_method,args.max_cat,comp1[i])

        # final report
        report(args.output_dir,args.project_name,args.pool_name)
    

# sentinel
if __name__ == "__main__":
    main()
