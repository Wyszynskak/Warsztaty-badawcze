library(shiny)
library(ggplot2)
library(survMisc)

cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG",
                "LUAD","LUSC","OV","READ","UCEC")
# Define UI for application that draws a histogram
shinyUI(fluidPage(
   img(src="t.png"),
   # Application title
   headerPanel(title="", windowTitle="Przglądarka mRNA"),
   #titlePanel("Aplikacja biomarkerów"),
   HTML('<font size="1">Autorzy: Karolina Wyszyńska, Krzysztof Rudaś</font>'),
   br(),
   br(),
   
   # Sidebar with a slider input for the number of bins
   sidebarLayout(
      sidebarPanel(
         selectInput("rak",
                     label  = h3(strong("Typ raka")),
                     choices  = cancer_types,
                     selected = cancer_types[1]
         ), 
         uiOutput("istotne_geny"),

         radioButtons("typ",
                      label  = h3(strong("Typ wykresu")),
                      choices  = c("Histogram", 
                                      "Krzywa przeżycia"),
                      selected = "Histogram"
         )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("wykres")
      )
   )
))