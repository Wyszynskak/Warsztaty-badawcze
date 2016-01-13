library(shiny)
library(ggplot2)
library(survMisc)
cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP",
                "LGG","LUAD","LUSC","OV","READ","UCEC")

for(i in 1:length(cancer_types)){
   eval(parse(text=paste0(cancer_types[i],"geny <- readRDS('../",cancer_types[i],"/geny.RDS')")))
   eval(parse(text=paste0(cancer_types[i],"histy <- readRDS('../",cancer_types[i],"/histy.RDS')")))
   eval(parse(text=paste0(cancer_types[i],"survy <- readRDS('../",cancer_types[i],"/survy.RDS')")))
}

shinyServer(function(input, output){
   output$istotne_geny<-renderUI({
      
      eval(parse(text=paste0("geny <- ", input$rak, "geny")))
      text <- paste('"', geny, '"', sep="")
      text <- paste(text, text, sep="=", collapse=", ")
      text <- paste0('list(', text, ')')
      lista <- eval(parse(text=text))
      selectInput("wybrany_gen", label=h3(strong("Istotne geny")), 
                  choices=lista, selected=lista[[1]])
      
   })
   output$wykres <- renderPlot({
      if(input$typ=="Histogram"){
         eval(parse(text=paste0("q <- ", input$rak, "histy[['",input$wybrany_gen,"']]")))
         q<-ggplot(q,aes_string(names(q)[1]))+ geom_histogram(binwidth = 0.1)+
            scale_y_continuous(name="Ilość oberwacji")+
            scale_x_continuous(name="Wartości")+
            ggtitle(paste0("Histogram mRNA ",input$wybrany_gen))+
            theme(plot.title = element_text(lineheight=20, face="bold"))
      } else {
         eval(parse(text=paste0("q <- ", input$rak, "survy[['",input$wybrany_gen,"']]")))
         eval(parse(text=paste0("l<-'krzywa przeżycia_dla_",input$wybrany_gen,"'")))
         q<-autoplot(q,xLab="czas",yLab="1-P(wystąpienie zdarzenia)",
                     title=paste0("Krzywa przeżycia dla mRNA ",input$wybrany_gen),
                     legTitle = "Warstwa",
                     legLabs=c("chorzy o wart. > mediana", "chorzy o wart. < mediana"))
      }
      q
   })
})