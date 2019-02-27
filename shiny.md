# Shiny app
<p align="justify"> To observe the results of the analysis in an interactive form, it is possible to launch a shiny app using the following command, with the same path used in the analysis: </p>

```R -e "shiny::runApp('/path/to/shiny/app/directory/',port=8080,host='localhost')" --args "/output_dir/project_name/pool_name"```

<p align="justify"> The shiny app will be available at the address indicated (in the example it is http://localhost:8080). It shows the results of the RNA-Seq analysis divided into a series of tabs for each phase. </p>

## Summary
<p align="center"><img src="/images/shiny1_summary.png" width="75%"></p>
<p align="justify"> The summary contains two tables. The first table shows the parameters chosen for the analysis, while the second shows a summary for the various samples. In particular, for each sample is reported the type (control or treated), the number of raw reads, which are the reads in the initial FastQ file, the number of Phix reads, which are the reads removed because they are part of Phix contaminating genome, and the number of ribosomal reads, which are instead the reads removed because belonging to the ribosomal DNA. </p>
<p align="justify"> In sidebar there is a button that allows to download all the results in a report in pdf, in html or in word. </p>

## FastQ quality
<p align="center"><img src="/images/shiny2_quality.png" width="75%"></p>
<p align="justify"> The FastQ quality tab contains the plots obtained from the quality analysis carried out with FastQ Screen and FastQC on the various samples. The first plot is the output of FastQ Screen and shows the percentage of DNA of sample reads mapped on different genomes and in particular human, murine, PhiX and ribosomal genomes. To calculate this percentage, 100,000 random reads are selected from the sample and they are sequenced on the reference genomes, then the result is multiplied by the rest of the genome. The second plot is instead the output of FastQC and shows the quality of the reads of the sample calculated using the Phred Score contained in the FastQ files. The drop-down menu in sidebar allows to browse through the various samples. </p>

## BAM quality
<p align="center"><img src="/images/shiny3_BAMquality.png" width="75%"></p>
<p align="justify"> The quality of the BAM files was evaluated using RSeQC software. In particular, in the first table are reported statistics related to the reads mapping, in the second table are reported reads fractions mapped on the coding exon part, on the 5'-UTR region, on the 3'-UTR region and on the intronic or intragenic regions, the final image shows the distribution of the internal distance between paired reads. </p>

## Differential expression analysis
<p align="center"><img src="/images/shiny4_DEA.png" width="100%"></p>
<p align="justify"> For the differential expression analysis, a summary table with the results and FPKM values for each sample is given. Then follows a series of summary plots and in particular: a PCA to evaluate the differences between the samples; a volcano plot, which reports the values of Fold Change and p-value for all genes; a heatmap of the 35 genes with greater variance, in which the value of the Z-score is reported and therefore the distance from the mean for the various samples; a heatmap showing the distances between the samples, calculated in a distance matrix using the euclidean distance method. </p>

## Meta-analysis
### GO analysis
<p align="center"><img src="/images/shiny5_GO.png" width="100%"></p>
<p align="justify"> The outputs of the Gene Ontology analysis are two summary tables, the first one showing a GO term for each line while the second one a gene for each line. Then follow three interactive networks, one for each category of GO, which allow to view the enriched genes and their Fold Change. Finally, there are a treemap, where the size of each rectangle is proportional to the number of genes, and three dotplots, which report the five terms of GO that were more enriched for the genes.</p>

### Pathway analysis
<p align="center"><img src="/images/shiny6_path.png" width="100%"></p>
<p align="justify"> For the pathway analysis a summary table with enriched pathways is shown, followed by an interactive network similar to that of GO. Finally, there is a dotplot with the five most enriched pathways. </p>
<p align="center"><img src="/images/shiny7_path1.png" width="100%"></p>
