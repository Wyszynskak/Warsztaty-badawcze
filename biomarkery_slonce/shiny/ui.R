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
# cancer_types<-c("READ","UCEC")
# cancer_names <- c("Rectum adenocarcinoma",
#                   "Uterine Corpus Endometrial Carcinoma")


cancer_list <- eval(parse(text=paste0("list(",paste('"',cancer_names,'"','=',
                                                    '"',cancer_types,'"',
                                              sep="", collapse=", "),")")))


shinyUI(fluidPage(
   img(src="english_title.png"),
   # Application title
   headerPanel(title="", windowTitle="mRNA Browser"),
   
   tabsetPanel(tabPanel("Instructions", 
                        helpText("hello world!"),
                        br(),
                        HTML('<font size="1">Version: 3.0 (18 January 2016)</font>'),
                        br(),
                        HTML('<font size="1">Authors: Karolina Wyszyńska, Krzysztof Rudaś</font>'),
                        br()),
   tabPanel("Application",
   # Sidebar with a slider input for the number of bins
   sidebarLayout(
      sidebarPanel(
         selectInput("cancer",
                     label  = h3(strong("Cancer type")),
                     choices  = cancer_list,
                     selected = cancer_list[1]
         ), 
         uiOutput("significant_genes"),

         radioButtons("typ",
                      label  = h3(strong("Plot type")),
                      choices  = c("Histogram", 
                                      "Survival curve"),
                      selected = "Histogram"
         ),
         conditionalPanel(condition="input.typ!='Histogram'",
                          uiOutput("plot_help_surv1"), br(),br())
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("wykres"),
         conditionalPanel(condition="input.typ=='Histogram'",
         uiOutput("plot_help_hist")),
         conditionalPanel(condition="input.typ!='Histogram'", br(),
         uiOutput("plot_help_surv2"),br(),br(),br())
      )),
   id="choice", selected="description"
   )
   )
))