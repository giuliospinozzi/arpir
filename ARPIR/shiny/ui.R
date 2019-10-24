library(shiny)
library(visNetwork)

args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

gen_sum=read.csv(paste0(out_dir,"/reports/general_report.csv"),sep="\t",header = F,row.names = 1)
sam_sum=read.csv(paste0(out_dir,"/reports/sample_report.csv"),sep="\t")
sam=as.character(sam_sum$sample_name)
comp=strsplit(as.character(gen_sum[7,1]),",")[[1]]

shinyUI(fluidPage(
  
  fluidRow( 
    shinyjs::useShinyjs(),
    singleton(tags$head(HTML(
      '
      <script type="text/javascript">
      $(document).ready(function() {
      // disable download at startup. download is the id of the downloadButton
      $("#download").attr("disabled", "true").attr("onclick", "return false;");
      
      Shiny.addCustomMessageHandler("download_ready", function(message) {
      $("#download").removeAttr("disabled").removeAttr("onclick").html(
      "<i class=\\"fa fa-download\\"></i>Generate report");
      });
      
      Shiny.addCustomMessageHandler("download_stop", function(message) {
      $("#download").attr("disabled", "true").attr("onclick", "return false;").html(
      "<i class=\\"fa fa-download\\"></i>Download");
      });
      
      })
      </script>
      '
    ))),
    navbarPage("RNA-Seq analysis",
               tabPanel("Summary", 
                        fluidRow(
                          column(2, 
                                 wellPanel(
                                   actionButton("start_proc", HTML("Click to start <br> processing data")),
                                   helpText("Download will be available once the processing is completed."),
                                   hr(),
                                   radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                                                inline = F),
                                   actionButton("download","Generate report"),
                                   br(),
                                   "It could require few minutes..."
                                 )
                          ),
                          column(9,align="center",
                                 textOutput("text"),
                                 tags$head(tags$style("#text{color: red;
                                                      font-size: 14px;
                                                      }"
                                 )
                                 ),
                                 br(),
                                 tags$b("Summary table"),
                                 br(),
                                 br(),
                                 tableOutput("gen_sum"),
                                 br(),
                                 tableOutput("sam_sum"))
                        )),
               tabPanel("FastQ quality",
                        br(),
                        fluidRow(
                          column(2,offset=0,
                                 wellPanel(
                                   selectInput("FastQ", "Choose a sample:", choices = sam)
                                 )
                          ),
                          column(8,align="center",
                                 plotOutput("quality1",height = "600px"),
                                 if (gen_sum[8,1]=="Paired_end") {
                                   plotOutput("quality2",height = "600px")
                                 },
                                 plotOutput("quality3",height = "600px"),
                                 if (gen_sum[8,1]=="Paired_end") {
                                   plotOutput("quality4",height = "600px")
                                 }
                          )
                        )
               ),
               tabPanel("BAM quality",
                        br(),
                        fluidRow(
                          column(2,offset=0,
                                 wellPanel(
                                   selectInput("BAM", "Choose a sample:", choices = sam)
                                 )
                          ),
                          column(8,
                                 verbatimTextOutput("stat"),
                                 verbatimTextOutput("read"),
                                 if (gen_sum[8,1]=="Paired_end") {
                                   plotOutput("inn",height = "800px")
                                 },
                                 plotOutput("jun_out",height = "800px"),
                                 plotOutput("spli_ev_out",height = "800px"),
                                 plotOutput("spli_jun_out",height = "800px")
                          )
                        )
               ),
               tabPanel("Differential expression analysis",
                        fluidRow(
                          column(6,offset=3,
                                 wellPanel(
                                   selectInput("comp", "Choose a comparison:", choices = comp)
                                 )
                          ),
                          br(),
                          column(12,offset=0,
                                 uiOutput("DEA_tab"),
                                 plotOutput("pca", height = "800px"),
                                 plotOutput("volcano", height = "800px"),
                                 plotOutput("heat1", height = "800px"),
                                 plotOutput("heat2", height = "800px")
                          )
                        )
               ),
               tabPanel("Meta-analysis",
                        fluidRow(
                          column(6,offset=3,
                                 wellPanel(
                                   selectInput("comp1", "Choose a comparison:", choices = comp)
                                 )
                          )
                        ),
                        uiOutput("Meta_analysis"),
                        if (!file.exists(paste0(out_dir,"/",strsplit(as.character(gen_sum[7,1]),",")[[1]][1],"/Meta-analysis"))){
                          column(12,align="center",
                                 br(),
                                 tags$span(style="font-size: 34px", strong("No meta-analysis"))
                          )
                        }
               ),
               tabPanel("References",
                        fluidRow(
                          tags$div(
                            tags$h1("References")
                          ),
                          tags$ul(style="font-size: 20px",
                                  tags$li("S. Andrews, \"FastQC: A quality control tool for high throughput sequence data.,\" Http://Www.Bioinformatics.Babraham.Ac.Uk/Projects/Fastqc/, 2010."), 
                                  tags$li("S. Andrews, \"FastQ Screen,\" http://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/, 2011."), 
                                  tags$li("H. Li and R. Durbin, \"Fast and accurate short read alignment with Burrows-Wheeler transform,\" Bioinformatics, vol. 25, no. 14, pp. 1754-1760, 2009."),
                                  tags$li("H. Li, \"A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data,\" Bioinformatics, vol. 27, no. 21, pp. 2987-2993, 2011."),
                                  tags$li("D. Kim, G. Pertea, C. Trapnell, H. Pimentel, R. Kelley, and S. L. Salzberg, \"TopHat2: Accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions,\" Genome Biol., vol. 14, no. 4, 2013."),
                                  tags$li("L. Wang, S. Wang, and W. Li, \"RSeQC: quality control of RNA-seq experiments,\" Bioinformatics, vol. 28, no. 16, pp. 2184-2185, 2012."),
                                  tags$li("D. Kim, B. Langmead, and S. L. Salzberg, \"HISAT: A fast spliced aligner with low memory requirements,\" Nat. Methods, vol. 12, no. 4, pp. 357-360, 2015."),
                                  tags$li("Y. Liao, G. K. Smyth, and W. Shi, \"FeatureCounts: An efficient general purpose program for assigning sequence reads to genomic features,\" Bioinformatics, vol. 30, no. 7, pp. 923-930, 2014."),
                                  tags$li("C. Trapnell, B. A. Williams, G. Pertea, A. Mortazavi, G. Kwan, M. J. Van Baren, S. L. Salzberg, B. J. Wold, and L. Pachter, \"Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation,\" Nat. Biotechnol., vol. 28, no. 5, pp. 511-515, 2010."),
                                  tags$li("L. A. Goff, C. Trapnell, and D. Kelley, \"cummeRbund: Analysis, exploration, manipulation, and visualization of Cufflinks high-throughput sequencing data,\" R Packag. Version 2.2, 2012."),
                                  tags$li("M. D. Robinson, D. J. McCarthy, and G. K. Smyth, \"edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.,\" Bioinformatics, vol. 26, no. 1, pp. 139-40, 2010."),
                                  tags$li("D. J. McCarthy, Y. Chen, and G. K. Smyth, \"Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation,\" Nucleic Acids Res., vol. 40, no. 10, pp. 4288-4297, 2012."),
                                  tags$li("M. I. Love, W. Huber, et al., \"Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2,\" Genome Biol., vol. 15, no. 12, p. 550, 2014."),
                                  tags$li("Yuan Tang, Masaaki Horikoshi, and Wenxuan Li. \"ggfortify: Unified Interface to Visualize Statistical Result of Popular R Packages.\" The R Journal 8.2, pp. 478-489, 2016."),
                                  tags$li("Masaaki Horikoshi and Yuan Tang, \"ggfortify: Data Visualization Tools for Statistical Analysis Results.\", 2016"),
                                  tags$li("J. Graffelman, \"calibrate: Calibration of Scatterplot and Biplot Axes.,\" R Packag. version 1.7.2., 2013."),
                                  tags$li("R. Gentleman, V. Carey, W. Huber, and F. Hahne, \"genefilter: methods for filtering genes from high-throughput experiments.,\" R Packag. version 1.60.0., 2017."),
                                  tags$li("G. R. Warnes, B. Bolker, L. Bonebakker, R. Gentleman, W. H. A. Liaw, T. Lumley, M. Maechler, A. Magnusson, S. Moeller, M. Schwartz, and B. Venables, \"gplots: Various R Programming Tools for Plotting Data,\" R Packag. version 3.0.1., 2016."),
                                  tags$li("Yihui Xie, \"DT: A Wrapper of the JavaScript Library 'DataTables'.\" R package version 0.4., 2018."),
                                  tags$li("Jeroen Ooms, \"magick: Advanced Graphics and Image-Processing in R.\" R package version 1.9., 2018."),
                                  tags$li("Almende B.V., Benoit Thieurmel and Titouan Robert, \"visNetwork: Network Visualization using 'vis.js' Library.\" R package version 2.0.4., 2018."),
                                  tags$li("Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson, \"shiny: Web Application Framework for R.\" R package version 1.0.5., 2017.")
                          )
                        ))
               
    ))
))

