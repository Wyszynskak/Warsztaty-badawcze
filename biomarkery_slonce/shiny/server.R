library(shiny)
library(ggplot2)
library(survMisc)
#All cancer types that we consider
cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP",
                "LGG","LUAD","LUSC","OV","READ","UCEC")
# cancer_types<-c("READ","UCEC")

#Read in data and pvalues for each combination (cancer,gene)
for(i in 1:length(cancer_types)){
   eval(parse(text=paste0(cancer_types[i],"data <- readRDS('./data/",
                          cancer_types[i],"dane.RDS')")))
#    eval(parse(text=paste0(cancer_types[i],"pvalue <- readRDS('./data/",
#                           cancer_types[i],"pvalue.RDS')")))
   eval(parse(text=paste0(cancer_types[i],"genes_list <- readRDS('./data/",
                          cancer_types[i],"genes_list.RDS')")))
}

shinyServer(function(input, output){
   #Creating select list with genes
   output$significant_genes<-renderUI({
      #Creating list to show as : "GENE_NAME (p-value: VALUE)"
      #Genes are already sorted from the most significat at the begining
      #If the value wasn't counted pvalue is "-"
      genes_list <- eval(parse(text=paste0("tmp <- ",input$cancer, 
                                           "genes_list")))
      selectInput("choosen_gene", label=h3(strong("Genes")), 
                  choices=tmp, selected=tmp[[1]])
      
   })
   #Creating stratum value slider for survival curve plot
   output$plot_help_surv1 <- renderUI({
     eval(parse(text=paste0("gene <- ",input$cancer, 
                            "data[, '",input$choosen_gene,"']")))
      max <- max(gene, na.rm=TRUE)
      min <- min(gene, na.rm=TRUE)
      med <- median(gene, na.rm=TRUE)
      sliderInput("stratum_value_choosen", label="Stratum value: ", min=min,
                    max=max,value=med, step=0.01, round=-2, ticks=FALSE)

   })
   
   #Creating text information under surv - plot
   output$plot_help_surv2 <- renderText({
     eval(parse(text=paste0("gene <- ",input$cancer, 
                            "data[, '",input$choosen_gene,"']")))
     n_low <- length(gene[gene<=input$stratum_value_choosen])
     n_high <- length(gene[gene>input$stratum_value_choosen])
      
     paste0("Number of observations <=", input$stratum_value_choosen,
                     ": ", n_low, br() , 
                     "Number of observations >", input$stratum_value_choosen,
                     ": ", n_high, br(),br(),
                     strong("Attention:"), br(),
                     "If the number of observations is small in one of 
                      above groups
                     the mRNA can seem to be significant, 
                     but it does not have to be.",br(),br())
     
   })
   
   #Create text information under histogram plot
   output$plot_help_hist <- renderUI({
     eval(parse(text=paste0("gene <- ",input$cancer, 
                            "data[, '",input$choosen_gene,"']")))
     helpText(paste0("Number of observations: ",length(gene[!is.na(gene)])))
     
   })
   
   #Creating plot itself
   output$wykres <- renderPlot({
      if(input$typ=="Histogram"){
        eval(parse(text=paste0("q <- ",
                               input$cancer, "data[,'",input$choosen_gene,
                               "',drop=FALSE]")))
         q<-ggplot(q,aes_string(names(q)[1]))+ geom_histogram(binwidth = 0.1)+
            scale_y_continuous(name="Number of observations")+
            scale_x_continuous(name="Values")+
            ggtitle(paste0("mRNA histogram ",input$choosen_gene))+
            theme(plot.title = element_text(lineheight=20, face="bold"))
      } else {
         eval(parse(text=paste0("tmp <- which(colnames(",input$cancer,
                                "data)=='", input$choosen_gene,"')")))
         eval(parse(text=paste0("dane <- ",input$cancer, "data[,c(1,",
                                tmp,",(ncol(",
                                input$cancer,"data)-1):ncol(",
                                input$cancer,"data))]")))
         dane$time <- dane$time/365
         med <- input$stratum_value_choosen
         eval(parse(text=paste0(
              'q<-survfit(Surv(time, status == "dead")~(', input$choosen_gene, 
                ">med), data=dane)")))
         try({suppressWarnings({
         q<-autoplot(q,xLab="years",yLab="",
                     title=paste0("Survival plot for mRNA ",input$choosen_gene,
                                  "\n1-P(death)"),
                     legTitle = "Stratum",
                     legLabs=c(paste0("Patients with value < ",round(med,4)),
                               paste0("Patients with value > ",
                                      round(med,4))))$plot+scale_y_continuous(
                                        limits=c(0,1))})}, silent=TRUE)
      }
      q
   })
})


