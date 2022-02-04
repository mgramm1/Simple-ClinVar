library(shiny)
library(ggplot2)
library(readr) 
library(stringr)
library(tidyr)
library(Hmisc)
library(ggrepel)
library(DT)
library(ggrepel)
library(dplyr)
library(tm)
library(shinycssloaders)
library(shinythemes)
library(drawProteins)
library(magrittr)
library(writexl)
library(plotly)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  navbarPage(theme = shinytheme("yeti"), id = "inTabset",selected = "panel1",
             title = "Simple ClinVar",
             ## Esta wea nose como funciona
             tags$head(tags$style(
               type="text/css",
               "#image1 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             ##
             tags$head(tags$style(
               type="text/css",
               "#image21 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             tags$head(tags$style(
               type="text/css",
               "#image22 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             tags$head(tags$style(
               type="text/css",
               "#image23 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             tags$head(tags$style(
               type="text/css",
               "#image31 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             tags$head(tags$style(
               type="text/css",
               "#image32 img {max-width: 100%;
               width: 100%;
               height: 100%;
               object-fit: contain;}"
             )),
             # Main box search and description -----------------------------------------
             tabPanel("Home", value = "panel1",
                      br(),br(),
                      h2("Simple ClinVar", align = "center"),
                      br(),
                      p("The content on this website is based on ",a("ClinVar", target="_blank", href="https://www.ncbi.nlm.nih.gov/clinvar/"), "database version July 14, 2021", align = "center"),
                      br(),
                      fluidRow(
                        column(width = 6, offset = 3, align = "center",
                               wellPanel("Type 'clinvar', a 'disease term', a 'gene name' or a 'variant' in HGVS format",
                                         style = "background-color: #333333;
                                         color: white;
                                         border-top-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-bottom: 0px;
                                         padding:5px"), 
                               wellPanel(br(),br(),
                                         textInput(inputId = "genename", label = NULL),#, value = "clinvar"),
                                         br(),
                                         actionButton(inputId ="goButton", label = "Submit", class = "btn-primary"), 
                                         style = "background-color: #ffffff;
                                         border-bottom-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-top: 0px")
                               ) # WellPanel
                               ), #Fluid row
                      fluidRow(column(width = 6, offset = 3, br(), br(), p("Simple ClinVar was developed to 
                                                                           provide gene- and disease-wise summary statistic based on all available genetic variants from ClinVar.\n
                                                                           How many missense variants are associated to heart disease?
                                                                           What are the top 10 genes mutated in Alzheimer? Does CDKL5 have pathogenic mutations? If so, where?
                                                                           Simple ClinVar is able to answer these questions and more, in a matter of seconds. \n ",
                                                                           align = "center"), p(""), style = "background-color: #ffffff")
                      ), br(),
                      fluidRow(column(width = 6, offset = 3, align = "center",imageOutput("imageNCBI")))
                      ),
             
             # Filtering "left" side -------------------------------------------------------------
             tabPanel("Results", value = "panel2",
                      fluidRow(column(width = 12, h6(strong("Displaying results for:")))), 
                      fluidRow(column(width = 2, style = "padding-right: 0px",
                                      wellPanel(textOutput(outputId = "displayquery"),
                                                style = "background-color:#008CBA;
                                                color:white;
                                                border-color:#008CBA;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-bottom: 0px;
                                                padding:15px;
                                                width: 80%"),
                                      fluidRow(column(width = 12,h5(""))),
                                      fluidRow(column(width = 12,h5(""))),
                                      wellPanel(h6(strong("Choose your filter")),
                                                style = "background-color:white;
                                                border-bottom: 2px solid #EEEEEE;
                                                border-top-color: white;
                                                border-right-color: white;
                                                border-left-color: white;
                                                box-shadow: 0px 0px 0px white;
                                                padding:3px;
                                                width: 100%"),
                                      actionButton(inputId = "filter", label = "Filter", width = "80%", 
                                                   style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                   margin-bottom: 0px;
                                                   padding:5px"),
                                      br(),br(),
                                     #tabsetPanel(type  = "tabs", id = "displayfilter", selected = "1", #display filters in tabs
                                     #              tabPanel(title = "1.", value = "1",
                                     #                       br(),
                                     #box(status = "primary", div(style = 'overflow-y: scroll; height:405px; width:200px; overflow: scroll',
                                                          checkboxGroupInput(inputId = "clinicalinput",
                                                                              label = "Clinical Significance",
                                                                              choices = c("Pathogenic",
                                                                                          "Likely pathogenic",
                                                                                          "Risk factor/Association",
                                                                                          "Uncertain/conflicting",
                                                                                          "Likely benign",
                                                                                          "Benign",
                                                                                          "Protective"),
                                                                              selected = NULL ),
                                                           checkboxGroupInput(inputId = "consequenceinput",
                                                                              label = "Consequence",
                                                                              choices= c("Missense",
                                                                                         "Stop gain",
                                                                                         "In frame indel",
                                                                                         "Frameshift",
                                                                                         "Synonymous",
                                                                                         "Splice-D/A",
                                                                                         "3-UTR",
                                                                                         "5-UTR",
                                                                                         "Non-coding",
                                                                                         "Intronic"),
                                                                              selected = NULL ),
                                                  #),
                                                  #tabPanel(title = "2.", value = "2",
                                                  #         br(),
                                                           checkboxGroupInput(inputId = "typeinput",
                                                                              label = "Type",
                                                                              choices= c("SNV",
                                                                                         "Indel",
                                                                                         "Insertion",
                                                                                         "Inversion",
                                                                                         "Deletion",
                                                                                         "Duplication",
                                                                                         "CNV loss",
                                                                                         "CNV gain",
                                                                                         "Complex"),
                                                                              selected = NULL ),
                                                           checkboxGroupInput(inputId = "reviewinput",
                                                                              label = "Review status",
                                                                              choices= c("Criteria provided/ multiple submitters/ no conflicts",
                                                                                         "Reviewed by expert panel",
                                                                                         "Practice guideline",
                                                                                         "Criteria provided/ single submitter",
                                                                                         "Criteria provided/ conflicting interpretations",
                                                                                         "No assertion criteria provided"),
                                                                              selected = NULL ),
                                                  #),
                                                  #tabPanel(title = "3.", value = "3",
                                                  #         br(),
                                                           checkboxGroupInput(inputId = "cadd",
                                                                              label = "CADD score",
                                                                              choices= c(">30",
                                                                                         ">25",
                                                                                         ">20"),
                                                                              selected = NULL ),
                                                           checkboxGroupInput(inputId = "gnomad",
                                                                              label = "Variant in gnomAD",
                                                                              choices= c("Yes",
                                                                                         "No"),
                                                                              selected = NULL )
                                                  #))
                                      #), #display filters in tabs
                                      # br(),br()
                                      ),
                               
                               # Results/Plots Panel -----------------------------------------------------
                               column(width = 10, style = "padding-left: 0px",
                                      fluidRow(
                                        column(2, offset = 1, 
                                               actionButton(inputId = "seevariants", label = textOutput(outputId = "displayobservation"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-success")),
                                        column(2, offset = 1,
                                               actionButton(inputId = "seegenes", label = textOutput(outputId = "displayngenes"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-danger")),
                                        column(2, offset = 1,
                                               actionButton(inputId = "seephenotypes", label = textOutput(outputId = "displayphenotypes"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-warning")),
                                        column(2, offset = 1,  
                                               actionButton(inputId = "seetable", label = "Show Table",
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px"))
                                               ),
                                      br(),br(),
                                      fluidRow(column(width = 12,
                                                      tabsetPanel(type  = "tabs", id = "displaywhat", selected = "1",
                                                                  tabPanel(title = NULL, value = "1",#br(),
                                                                           fluidRow(
                                                                             column(3, h5(strong("\tType")), align = "center",
                                                                                    withSpinner(plotOutput("Type", width="100%"))),
                                                                             column(3, h5(strong("\tConsequence")), align = "center",
                                                                                    withSpinner(plotOutput("consequence", width="100%"))),
                                                                             column(3, h5(strong("\tClinical Significance")), align = "center",
                                                                                    withSpinner(plotOutput("ClinicalSignificance", width="100%"))),
                                                                             column(3, h5(strong("\tReview status")), align = "center",
                                                                                    withSpinner(plotOutput("review", width="100%"))))
                                                                  ), 
                                                                  tabPanel(title = NULL, value = "2",
                                                                           fluidRow(
                                                                           uiOutput("genefield")
                                                                           )
                                                                  ),
                                                                  tabPanel(title = NULL, value = "3",#br(),
                                                                           fluidRow(column(8, h5(strong("Top 10 phenotypes associated")), align = "center", br(),
                                                                                           plotOutput("DiseaseCount")),
                                                                                    column(4, h5(strong("Total phenotypes associations")), align = "center",
                                                                                           wellPanel(DT::dataTableOutput("DiseaseCountTotal"),
                                                                                                     style = "background-color: #ffffff;
                                                                                                     border-color: #ffffff;
                                                                                                     box-shadow: 0px 0px 0px #ffffff;
                                                                                                     margin-bottom: 5px")
                                                                                           )
                                                                                           ) #Fluid row
                                                                  ),
                                                                  tabPanel(title = NULL, value = "4",
                                                                           fluidRow(column(12,br(), align = "center",
                                                                                           wellPanel(DT::dataTableOutput("clinvartable1"),
                                                                                                     style = "background-color: #ffffff;
                                                                                                     border-color: #ffffff;
                                                                                                     box-shadow: 0px 0px 0px #ffffff;
                                                                                                     margin-bottom: 5px")
                                                                                           ) # Column
                                                                                           ), 
                                                                           fluidRow(column(12,downloadLink("downloadData1", "Download (.csv)"), downloadLink("downloadData1_xlsx", "Download (.xlsx)")
                                                                                           )), 
                                                                           fluidRow(column(12,br(), align = "center",
                                                                                          wellPanel(DT::dataTableOutput("clinvartable2"),
                                                                                                   style = "background-color: #ffffff;
                                                                                                   border-color: #ffffff;
                                                                                                   box-shadow: 0px 0px 0px #ffffff;
                                                                                                   margin-bottom: 5px")
                                                                                                              )# Column
                                                                                                              ),# Fluid row
                                                                           fluidRow(column(12,downloadLink("downloadData2", "Download (.csv)"), downloadLink("downloadData2_xlsx", "Download (.xlsx)"))
                                                                           )
                                                                  ) #TabPanle
                                                      ) #TabsetPanel
                                                  ) # Column
                                              ), #Fluidrow
                                              fluidRow(column(width = 12,
                                                       wellPanel(style = "background-color:white;
                                                                          border-bottom: 2px solid #EEEEEE;
                                                                          border-top-color: white;
                                                                          border-right-color: white;
                                                                          border-left-color: white;
                                                                          box-shadow: 0px 0px 0px white;
                                                                          padding:0px;
                                                                          width: 100%"))), 
                                              fluidRow(column(width = 2, offset = 10,
                                                              actionButton(inputId = "newsearch", label = "New Search", 
                                                                            width = "80%",
                                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                                            margin-bottom: 0px;
                                                                            padding:5px")
                                                              )
                                              ) #Fluidrow 
                        ) # column de Results Panel
                      ),
                      fluidRow(column(width = 12,
                                      wellPanel(style = "background-color:white;
                                                border-bottom: 2px solid #EEEEEE;
                                                border-top-color: white;
                                                border-right-color: white;
                                                border-left-color: white;
                                                box-shadow: 0px 0px 0px white;
                                                padding:0px;
                                                 width: 100%"))) 
                    #,verbatimTextOutput("test")
             ), #mainTabPanel,
             
             # about page --------------------------------------------------------------
             tabPanel("About", value = "panel3",                  
                      br(),br(),
                      h2("Simple ClinVar", align = "center"),
                      br(),
                      p("For detailed information about", a("Simple ClinVar", target="_blank", href="http://simple-clinvar.broadinstitute.org"), "please refer to the original publication:", align = "center"),
                      br(),
                      fluidRow(column(width = 6, offset = 3, align = "center",
                                      wellPanel("Citation and contact",
                                                style = "background-color: #333333;
                                                color: white;
                                                border-top-color: #333333;
                                                border-left-color: #333333;
                                                border-right-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-bottom: 0px;
                                                padding:5px"), 
                                      wellPanel(br(),p("Eduardo Perez-Palma, Marie Gramm, Peter Nürnberg, Patrick May and Dennis Lal.",
                                                       a("Simple ClinVar: an interactive web server to explore and retrieve gene and disease variants aggregated in ClinVar database", target="_blank", href="https://doi.org/10.1093/nar/gkz411"),
                                                       strong(em(". Nucleic Acids Research")), "(2019) PMID:31114901."),
                                                br(),
                                                p(icon("envelope", lib = "glyphicon"),"   eduardoperez@udd.cl |",
                                                  icon("envelope", lib = "glyphicon"),"   mgramm1@smail.uni-koeln.de |",
                                                  icon("envelope", lib = "glyphicon"),"   lald@ccf.org "),
                                                style = "background-color: #ffffff;
                                                border-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-top: 0px")
                      ) # Column
                      ), #Fluid row
                      fluidRow(column(width = 6, offset = 3, align = "center",
                                      h3("Summary", align = "left"), 
                                      p("Clinical genetic testing has exponentially expanded in recent years, leading to an overwhelming amount of patient variants with high variability in pathogenicity and heterogeneous phenotypes. A large part of the variant level data are comprehensively aggregated in public databases such as ClinVar and are publicly accessible.  However, the ability to explore this rich resource and answer general questions such as “How many missense variants are associated to a specific disease or gene?” or “In which part of the protein are patient variants located?” is limited and requires advanced bioinformatics processing.", align = "justify"),
                                      p("Here, we present ", 
                                        a("Simple ClinVar", target="_blank", href="http://scv.broadinstitute.org/"), "a web-server application that is able to provide variant, gene, and disease level summary statistics based on the entire ClinVar database in a dynamic and user-friendly web-interface.  Overall, our web application is able to interactively answer basic questions regarding genetic variation and its known relationships to disease. Our website will follow ClinVar monthly releases and provide easy access to the rich ClinVar resource to a broader audience including basic and clinical scientists.",  align = "justify"),
                                      h3("Methods", align = "left"), 
                                      p("The ClinVar database is downloaded on a monthly basis directly from the", 
                                        a("ClinVar ftp site", target="_blank", href="ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/"),
                                        "The tabular data is processed internally to produce a pre-filtered ClinVar file for the user to explore on the Simple ClinVar web-server. The pre-filtering step is designed to reduce the complexity of ClinVar entries as well as to provide fast access to high quality entries.", align = "justify"),
                                      p("Detailed description of the pre-filtering step as well as access to the complete pre-filtering pipeline is provided at our",
                                        a("GitHub repository.",  target="_blank", href="http://github.com/dlal-group/scv"), 
                                      "Briefly the pipeline performs four main tasks:", align = "justify"),
                                      p(strong("1)	"),"First, we keep only entries from the human reference genome version GRCh37.p13/hg19 and referring to canonical transcripts.", align = "justify"),
                                      p(strong("2)	"),"Second, Molecular consequence is inferred through the analysis of the Human Genome Variation Society (HGVS) sequence variant nomenclature field. HGVS format contains string patterns that allow molecular consequence inference (e.g. ENST00001234 c.3281 [p.Gly1046Arg]). Specifically, when the variant is reported to cause an amino acid change different to the reference it is annotated as a “missense” variant (e.g. p.Gly1046Arg). If the genetic variant leads to the same amino acid (e.g. p.Gly1046Gly) or a stop codon (e.g. p.Gly1046*) the entry is annotated as “synonymous” or “stop-gain” variant, respectively. Depending on the observed outcome, small insertions and deletions molecular consequence (collectively “indels”; e.g. p.Gly1046Glyfs., p.Gly1046Glyins., p.Gly1046del.*) are separated in to “frameshifts” and “in-frame indels”.", align = "justify"),
                                      p(strong("3)	"),"Third, we reduced the complexity of the clinical significance field by regrouping and merging them in to five unique and non-redundant categories: “Pathogenic”, “Likely pathogenic”, “Risk factor and Association”, “Protective/Likely benign” and “Benign”. Conflicting interpretations of pathogenicity, variants of unknown significance (VUS) and contradictory evidence (e.g. “Likely benign” alongside “Risk” evidences) were combined together in to an “Uncertain/Conflicting” category. Similarly, variants annotated with multiple evidence categories of the same evidence direction such as “Pathogenic” alongside “Likely pathogenic” were combined and the respective lower evidence category assigned. Fourth, ClinVar entries with phenotypes annotated as “not provided” and “not specified” were combined in to one single category called “Not provided / Not specified”.", align = "justify"),
                                      p(strong("4)	"),"Fourth, ClinVar entries with missing annotations such as the absence of anHGVS variant name or incomplete genomic coordinates are filtered out. Currently, 493,240 out of 503,065 (98.04%) ClinVar entries (April 22 release) are included in Simple ClinVar.", align = "justify"),
                                      
                                      p("Interactive summary statistics, variant mapping and visualization was developed with the Shiny framework of R studio software. App deployment, hosting and update is performed with Google Cloud services (Figure 1).", align = "justify")
                                      )),
                      br(),
                      fluidRow(column(width = 6, offset = 3, align = "center",
                                      imageOutput("image1"),
                                      br(),br(),
                                      h3("Simple ClinVar main features", align = "left"),
                                      p("From the front page of Simple ClinVar the user can submit three types of queries: ", align = "justify"),
                                      p(strong("1)	Database-wise query:"),"Triggered by submitting without a query or with the keyword “clinvar”, it will yield summary statistics of the entire ClinVar database. By the time of submission (ClinVar February 2019 release) Simple ClinVar contains 493,240 genetic variants, identified in 18,502 genes found in patients with 11,098 phenotypes. The database query mode coupled with dynamic filtering allows the user to explore which are the most common disorders and types of variants most commonly found in the whole database. Similarly, it is possible to evaluate immediately which are the genes with the most pathogenic variations or variants of unknown significance (Figure 2).", align = "justify"),
                                      wellPanel(imageOutput("image21"), style= "border-color: #ffffff;background-color: #ffffff;"),
                                      p(strong("2) Gene-wise query:"),"Submitting a RefSeq gene name on the main page will forward the user to the corresponding summary statistic page for all the genetic variants annotated in the gene. For example, querying for",em("“CDKL5”")," will show 675 genetic variants currently associated with 27 genetic disorders. Here, the user will see these variants mapped over the corresponding protein sequence alongside domain information from UniProt. Furthermore, the user can explore unique gene-specific variant statistics such as determining how many clinical phenotypes are associated with a given gene or where the pathogenic versus benign variants are located in the protein (Figure 3). ", align = "justify"),
                                      wellPanel(imageOutput("image22"), style= "border-color: #ffffff;background-color: #ffffff;"),
                                      p(strong("3)	Disease-term-wise query:"),"Querying a broad disease term of interest in Simple ClinVar will provide the user with all genes, variants and phenotypes associated to the given disease term. As an example, the disease term query “heart” will yield 814 genetic variants in 61 genes associated with 233 phenotypes, with missense SNVs (n=321) as the most common variant type. Here, for each selected disease term the user can answer general questions such as: how many genes are associated with heart disease? How many annotated terms and disorder subtypes can be found related to heart disease? (Figure 4). ", align = "justify"),
                                      wellPanel(imageOutput("image23"), style= "border-color: #ffffff;background-color: #ffffff;"),
                                      p("Independently of the input mode, the output displayed at the results tab can always be explored between four sections marked by the top square buttons in the colors green, red, orange and grey. The user can switch between the color areas. The green button will show all genetic variants available for the query and see the counts of variant type, molecular consequence, clinical significance, and review status. The red and orange buttons will show the top ten genes and phenotypes associated with the corresponding query and the complete list in table mode, respectively. In the case of a gene-wise query submission, the red button will show the variant mapping over the canonical protein sequence of the gene queried. Finally, the grey button will show the table mode were the subset of the pre-filtered ClinVar file currently in display is shown and available for download for downstream analysis.", align = "justify"),
                                      br(),br(),
                                      h3("Simple ClinVar filtering examples", align = "left"),
                                      p("At all query levels, the output can be dynamically filtered by variant type, molecular consequence, clinical significance and review status in any combination. We show two examples of how this feature can be used as a fast and powerful tool for clinical researchers. ", align = "justify"),
                                      p(strong("Example 1:"),"We use “epilepsy” as a disease term query and display the red button gene view area. Unfiltered results show the top ten “epilepsy” genes associated in descending order according to the number of qualifying variants:",em("SCN1A, SCN9A, CACNA1H, GRIN2A, DEPDC5, RELN, KCNT1, KCNQ3, ALDH7A1")," and ", em("CHRNA4"),". Next, after filtering for pathogenic variants, the top ten gene list is updated to",em("SCN1A, DEPDC5, GRIN2A, SCARB2, ALDH7A1, LGI1, MEF2C, NPRL3, SCN9A,")," and ", em("SPATA5.")," The user can conclude that the order of frequently mutated genes and genes with the most pathogenic classified genes is not the same. In the example, only ",em("SCN1A, DEPDC5,")," and ",em("GRIN2A")," are in the top ten “epilepsy” gene list both as genes with the most variants and most pathogenic variants (Figure 5). ", align = "justify"),
                                      wellPanel(imageOutput("image31"), style= "border-color: #ffffff;background-color: #ffffff;"),
                                      p(strong("Example 2:"),"We evaluate the gene-wise query for ",em("“SCN2A”")," on the red button gene view. Currently, there are variants mapped on the protein sequence of",em("SCN2A."), "If we filter for “Missense” and “Pathogenic” only 40 variants remain and are concentrated inside the transmembrane domains. The user can conclude that these regions containing the majority of pathogenic variants are of key importance for the protein function (Figure 6).", align = "justify"),
                                      wellPanel(imageOutput("image32"), style= "border-color: #ffffff;background-color: #ffffff;")
                      )),
                      fluidRow(column(width = 12, br()))

                      ) # TabPanel
              ) # NavbarPage
            ))
