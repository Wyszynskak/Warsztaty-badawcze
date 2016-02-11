library(shiny)
library(ggplot2)
library(survMisc)

cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG",
                "LUAD","LUSC","OV","READ","UCEC")
cancer_names <- c("Breast invasive carcinoma", "Colon adenocarcinoma",
                  "Colon adenocarcinoma & Rectum adenocarcinoma",
                  "Glioblastoma multiforme & Brain Lower Grade Glioma",
                  "Kidney: Chromophobe, renal clear cell carcinoma, renal papillary cell carcinoma",
                  "Kidney renal clear cell carcinoma",
                  "Kidney renal papillary cell carcinoma","Brain Lower Grade Glioma",
                  "Lung adenocarcinoma","Lung squamous cell carcinoma",
                  "Ovarian serous cystadenocarcinoma","Rectum adenocarcinoma",
                  "Uterine Corpus Endometrial Carcinoma")

cancer_list <- eval(parse(text=paste0("list(",paste('"',cancer_names,'"','=',
                                                    '"',cancer_types,'"',
                                              sep="", collapse=", "),")")))


shinyUI(fluidPage(
   #Title an window title
   HTML('<div align="center">'),
   img(src="english_title.png"),
   HTML ('</div>'),
   headerPanel(title="", windowTitle="mRNA Browser"),
   
   #Instruction panel
   conditionalPanel(condition="input.instructionButton>=input.appButton",
     h2("Welcome in mRNA Browser!"),
     br(),
     "It is an application that helps exploring the data from",
     a("RTCA", href='https://tcga-data.nci.nih.gov/tcga/'),
     " project.", br(), "All the data was analased using ", 
     a("RTCGA package", href= "https://github.com/RTCGA"),
     "dedicated to this project.", br(),
     "For any further information see", a("documentation. ", 
       href="http://mi2.mini.pw.edu.pl:8080/RTCGA/KRWK/shiny/mRNABrowser-documentation.pdf"),
     br(),
     br(),
     HTML('<font size="1">Version: 4.0 (12 February 2016)</font>'),
     br(),
     HTML('<font size="1">Authors: Karolina Wyszyńska, Krzysztof Rudaś</font>'),
     br(),
     HTML('<div align="center">'),
     actionButton("appButton", "Go to Application!"),
     HTML('</div>')
   ),
   
   #Application panel
   conditionalPanel(condition="input.instructionButton<input.appButton",
      HTML('<div align="right">'),
      actionButton("instructionButton", "Instruction"),
      HTML('</div>'),
      
      tabsetPanel(
        tabPanel("Gene influence on diffrent types of cancer", value="compare",
                 fluidPage(
                   fluidRow(column(4),
                            column(8, 
                                   uiOutput("compareCancerTitle"))),
                   fluidRow(
                     column(4, class="well",
                            br(),

                            column(2, h5(strong("Gene: "))),
                            column(10,uiOutput("significant_genesCompGenes")),
                            HTML('<font size="1">Choose a gene which influence
                                  you wish to compare between different types
                                 of cancer. Remember that comparing is based
                                 on mRNA expression. </font>'), br(), br(),
                            br(),
                            h5(strong("Cancer types: ")),
                            fluidRow(
                              column(6,
                                     selectInput("cancerComp1",label  =NULL,
                                                 choices  = cancer_list,
                                                 selected = cancer_list[1])),
                              column(6,
                                     selectInput("cancerComp2",label  = NULL,
                                                 choices  = cancer_list,
                                                 selected = cancer_list[2]))
                            ),
                            fluidRow(
                              column(6,
                                     selectInput("cancerComp3",label  = NULL,
                                                 choices  = cancer_list,
                                                 selected = cancer_list[3])),
                              column(6,
                                     selectInput("cancerComp4",label  =NULL,
                                                 choices  = cancer_list,
                                                 selected = cancer_list[4]))
                            ),
                            
                            HTML('<font size="1">Choose cancer types 
                                 which you want to include in the comparison. </font>'),
                            br(),br(),br(),
                            fluidRow( column(2, h5(strong("Years: "))),
                                      column(10,  uiOutput("yearsValue"))
                            ),
                            HTML('<font size="1">Manipulating years axis for
                                 better comparison.</font>'),br()
                     ),
                     column(4, 
                            plotOutput("survplotComp1", width="100%"),
                            
                            plotOutput("survplotComp3", width="100%")
                            
                     ),
                     column(4,
                            plotOutput("survplotComp2", width="100%"),
                            
                            plotOutput("survplotComp4", width="100%")
                            
                     )
                   )
                 )
                 
        ),
        tabPanel("Genes influences on specyfic type of cancer", value="compareGene",
                 fluidPage(
                   fluidRow(column(4),
                            column(8, 
                                   uiOutput("compareGenesTitle"))),
                   fluidRow(
                     column(4, class="well",
                            br(),

                            column(2, h5(strong("Cancer: "))),
                            column(10,selectInput("choosen_CompCancer",
                                                  label=NULL, 
                                                  choices=cancer_list,
                                                  selected=cancer_list[[1]])
                            ),
                            HTML('<font size="1">Choose a cancer type on which influence
                                 you wish to investigate. </font>'),
                            br(),br(),br(),

                            h5(strong("Genes: ")),
                            fluidRow(
                              column(6,uiOutput("gene1")),
                              column(6,uiOutput("gene3"))
                            ),
                            fluidRow(
                              column(6,uiOutput("gene2")),
                              column(6,uiOutput("gene4"))
                            ),
                            HTML('<font size="1">Choose genes which influences
                                 you want to compare. Remember that results are based
                                 on mRNA expression. </font>'),
                            br()
                            ),
                     column(4, 
                            plotOutput("survplotGene1", width="100%"),
                            
                            plotOutput("survplotGene3", width="100%")
                            
                     ),
                     column(4,
                            plotOutput("survplotGene2", width="100%"),
                            
                            plotOutput("survplotGene4", width="100%")
                            
                     )
                     )
                   )
                 
        ),
          tabPanel("Detailed Survival Plot", value="SurvPlot",
                   # Sidebar with a slider input for the number of bins
                   sidebarLayout(
                     sidebarPanel(width=4,
                       br(),
                       fluidRow(
                         column(2,
                                h5(strong("Cancer: "))),
 
                         column(10,
                                selectInput("cancerSurv",label  =NULL,
                                            choices  = cancer_list,
                                            selected = cancer_list[1])
                         ), 
                         HTML('<font size="1">\tType of a cancer that patiens 
                              included on surivival plot were
                              suffering from.</font>')
                       ), br(),br(),
                       fluidRow(
                         column(2, 
                                h5(strong("Gene: "))),
                         column(10,
                                uiOutput("significant_genesSurv"))
                       ),
                       HTML('<font size="1">Gene influence measured by mRNA 
                                          expression. P-value was calculated
                            when median was taken for a stratum value.</font>'),
                       br(), br(),
                       uiOutput("plot_help_surv1"),
                       plotOutput("histogram", height="200px"),
                       br(),
                       tableOutput("adres"),
                       br(),
                       actionButton("save","Save Kaplan-Meier Plot") 
                     ),
                                                              
   
                     
                     # Show a plot of the generated distribution
                     mainPanel(width=8,
                        plotOutput("survplot"),
                        uiOutput("plot_help_surv2")
                     )
                   )
        ),
      id="typ", selected="compare" 
      )
  )
))