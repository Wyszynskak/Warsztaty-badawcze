library(shiny)
library(ggplot2)
library(survMisc)
#All cancer types that we consider
cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP",
                "LGG","LUAD","LUSC","OV","READ","UCEC")
cancer_names <- c("Breast invasive carcinoma", "Colon adenocarcinoma",
                  "Colon adenocarcinoma & Rectum adenocarcinoma",
                  "Glioblastoma multiforme & Brain Lower Grade Glioma",
                  "Kidney: Chromophobe, renal clear cell carcinoma, renal papillary cell carcinoma",
                  "Kidney renal clear cell carcinoma",
                  "Kidney renal papillary cell carcinoma","Brain Lower Grade Glioma",
                  "Lung adenocarcinoma","Lung squamous cell carcinoma",
                  "Ovarian serous cystadenocarcinoma","Rectum adenocarcinoma",
                  "Uterine Corpus Endometrial Carcinoma")

cancer_list1 <- eval(parse(text=paste0("list(",paste('"',cancer_types,'"','=',
                                                    '"',cancer_names,'"',
                                                    sep="", collapse=", "),")")))

#Read in data and pvalues for each combination (cancer,gene)
for(i in 1:length(cancer_types)){
   eval(parse(text=paste0(cancer_types[i],"data <- readRDS('./data/",
                          cancer_types[i],"dane.RDS')")))
   eval(parse(text=paste0(cancer_types[i],"genes_list <- readRDS('./data/",
                          cancer_types[i],"genes_list.RDS')")))
}


shinyServer(function(input, output){
    ############################################################################
    ############    Gene influence on diffrent types of cancer    ##############
    ############################################################################
  output$compareCancerTitle <- renderUI({
    h3(strong(paste0("Survival plots for ", input$choosen_geneCompGenes)),
       align="center")
  })
  
  output$significant_genesCompGenes<-renderUI({
    #Creating list to show as : "GENE_NAME (p-value: VALUE)"
    #Genes are already sorted from the most significat at the begining
    #If the value wasn't counted pvalue is "-"
     eval(parse(text=paste0("all_genes_list <- readRDS('./data/",
                            "ALLCancerGenes_order.RDS')")))
    selectInput("choosen_geneCompGenes", label=NULL, 
                choices=all_genes_list, selected=all_genes_list[[1]])
    
  })
  
  #Moving x axis values
  output$yearsValue <- renderUI({
    if(is.null(input$cancerComp1) || is.null(input$cancerComp2) ||
       is.null(input$cancerComp3) || is.null(input$cancerComp4) || 
       is.null(input$choosen_geneCompGenes)){
      return(    sliderInput("yearsValueSlider", label=NULL, 0, 20, 
                             value=c(0,20), step=0.5))
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp1,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp1, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp1,"data)-1):ncol(",
                           input$cancerComp1,"data))]")))
    dane$time <- dane$time/365
    
    maximum <- max(dane$time, na.rm=TRUE)
    
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp2,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp2, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp2,"data)-1):ncol(",
                           input$cancerComp2,"data))]")))
    dane$time <- dane$time/365
    
    maximum <- c(maximum, max(dane$time, na.rm=TRUE))
    
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp3,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp3, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp3,"data)-1):ncol(",
                           input$cancerComp3,"data))]")))
    dane$time <- dane$time/365
    
    maximum <- c(maximum, max(dane$time, na.rm=TRUE))
    
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp4,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp4, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp4,"data)-1):ncol(",
                           input$cancerComp4,"data))]")))
    dane$time <- dane$time/365
    
    maximum <- ceiling(max(maximum, max(dane$time, na.rm=TRUE), na.rm=TRUE))

    
    sliderInput("yearsValueSlider", label=NULL, 0, maximum, 
                value=c(0,maximum), step=0.5)
    
  })
  
  output$survplotComp1 <- renderPlot({
    if(is.null(input$cancerComp1) || is.null(input$choosen_geneCompGenes)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp1,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp1, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp1,"data)-1):ncol(",
                           input$cancerComp1,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_geneCompGenes, 
      ">med), data=dane)")))

    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(cancer_list1[[input$cancerComp1]],
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                     limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                             limits=input$yearsValueSlider, name="years"
                                           )
    })
    q
  })
  
  output$survplotComp2 <- renderPlot({
    if(is.null(input$cancerComp2) || is.null(input$choosen_geneCompGenes)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp2,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp2, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp2,"data)-1):ncol(",
                           input$cancerComp2,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_geneCompGenes, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(cancer_list1[[input$cancerComp2]],
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                     limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                             limits=input$yearsValueSlider, name="years"
                                           )
    })
    q
  })
  
  output$survplotComp3 <- renderPlot({
    if(is.null(input$cancerComp3) || is.null(input$choosen_geneCompGenes)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp3,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp3, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp3,"data)-1):ncol(",
                           input$cancerComp3,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_geneCompGenes, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(cancer_list1[[input$cancerComp3]],
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                     limits=c(0,1))+theme(
                                       legend.justification=c(0,0),
                                      legend.position=c(0,0))+scale_x_continuous(
                                        limits=input$yearsValueSlider, name="years"
                                      )
      })
    q
  })
  
  output$survplotComp4 <- renderPlot({
    if(is.null(input$cancerComp4) || is.null(input$choosen_geneCompGenes)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$cancerComp4,
                           "data)=='", input$choosen_geneCompGenes,"')")))
    eval(parse(text=paste0("dane <- ",input$cancerComp4, "data[,c(1,",
                           tmp,",(ncol(",
                           input$cancerComp4,"data)-1):ncol(",
                           input$cancerComp4,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_geneCompGenes, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(cancer_list1[[input$cancerComp4]],
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                     limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                             limits=input$yearsValueSlider, name="years"
                                           )
      })
    q
  })
  
  
  ############################################################################
  ############    Genes influences on specyfic type of cancer    ##############
  ############################################################################
  output$compareGenesTitle <- renderUI({
    h3(strong(paste0("Survival plots for ", 
                     cancer_list1[[input$choosen_CompCancer]])),
                      align="center")
  })
  
  output$gene1 <- renderUI({
    eval(parse(text=paste0("tmp <- ",input$choosen_CompCancer, 
                           "genes_list")))
    selectInput("choosen_gene1", label=NULL, 
                choices=tmp, selected=tmp[[1]])
  })
  output$gene2 <- renderUI({
    eval(parse(text=paste0("tmp <- ",input$choosen_CompCancer, 
                           "genes_list")))
    selectInput("choosen_gene2", label=NULL, 
                choices=tmp, selected=tmp[[2]])
  })
  output$gene3 <- renderUI({
    eval(parse(text=paste0("tmp <- ",input$choosen_CompCancer, 
                           "genes_list")))
    selectInput("choosen_gene3", label=NULL, 
                choices=tmp, selected=tmp[[3]])
  })
  output$gene4 <- renderUI({
    eval(parse(text=paste0("tmp <- ",input$choosen_CompCancer, 
                           "genes_list")))
    selectInput("choosen_gene4", label=NULL, 
                choices=tmp, selected=tmp[[4]])
  })

  output$survplotGene1 <- renderPlot({
    if(is.null(input$choosen_CompCancer) || is.null(input$choosen_gene1)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$choosen_CompCancer,
                           "data)=='", input$choosen_gene1,"')")))
    eval(parse(text=paste0("dane <- ",input$choosen_CompCancer, "data[,c(1,",
                           tmp,",(ncol(",
                           input$choosen_CompCancer,"data)-1):ncol(",
                           input$choosen_CompCancer,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_gene1, 
      ">med), data=dane)")))
    
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(input$choosen_gene1,
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                                                           limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                             name="years"
                                           )
    })
    q
  })
  
  output$survplotGene2 <- renderPlot({
    if(is.null(input$choosen_CompCancer) || is.null(input$choosen_gene2)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$choosen_CompCancer,
                           "data)=='", input$choosen_gene2,"')")))
    eval(parse(text=paste0("dane <- ",input$choosen_CompCancer, "data[,c(1,",
                           tmp,",(ncol(",
                           input$choosen_CompCancer,"data)-1):ncol(",
                           input$choosen_CompCancer,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_gene2, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(input$choosen_gene2,
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                                                           limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                              name="years"
                                           )
    })
    q
  })
  
  output$survplotGene3 <- renderPlot({
    if(is.null(input$choosen_CompCancer) || is.null(input$choosen_gene3)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
    eval(parse(text=paste0("tmp <- which(colnames(",input$choosen_CompCancer,
                           "data)=='", input$choosen_gene3,"')")))
    eval(parse(text=paste0("dane <- ",input$choosen_CompCancer, "data[,c(1,",
                           tmp,",(ncol(",
                           input$choosen_CompCancer,"data)-1):ncol(",
                           input$choosen_CompCancer,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_gene3, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(input$choosen_gene3,
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(
                                     name="",limits=c(0,1))+theme(
                                       legend.justification=c(0,0),
                                       legend.position=c(0,0))+scale_x_continuous(
                                       name="years"
                                       )
    })
    q
  })
  
  output$survplotGene4 <- renderPlot({    
    if(is.null(input$choosen_CompCancer) || is.null(input$choosen_gene4)){
    p <-plot.new()
    text(0.5,0.5,"LOADING...")
    return (p)
  }
    eval(parse(text=paste0("tmp <- which(colnames(",input$choosen_CompCancer,
                           "data)=='", input$choosen_gene4,"')")))
    eval(parse(text=paste0("dane <- ",input$choosen_CompCancer, "data[,c(1,",
                           tmp,",(ncol(",
                           input$choosen_CompCancer,"data)-1):ncol(",
                           input$choosen_CompCancer,"data))]")))
    dane$time <- dane$time/365
    med <- median(dane[ ,2], na.rm=TRUE)
    eval(parse(text=paste0(
      'q<-survfit(Surv(time, status == "dead")~(', input$choosen_gene4, 
      ">med), data=dane)")))
    suppressMessages({
      q<-autoplot(q,xLab="years",yLab="",
                  title=paste0(input$choosen_gene4,
                               "\n1-P(death)"),
                  legTitle = "Stratum",
                  legLabs=c(paste0("Patients with value < ",round(med,4)),
                            paste0("Patients with value > ",
                                   round(med,4))))$plot+scale_y_continuous(name="",
                                                                           limits=c(0,1)
                                   )+theme(legend.justification=c(0,0),
                                           legend.position=c(0,0))+scale_x_continuous(
                                             name="years"
                                           )
    })
    q
  })
  
  
    ############################################################################
    ############               Detailed survival plot             ##############
    ############################################################################  
  output$significant_genesSurv<-renderUI({
    #Creating list to show as : "GENE_NAME (p-value: VALUE)"
    #Genes are already sorted from the most significat at the begining
    #If the value wasn't counted pvalue is "-"
    eval(parse(text=paste0("tmp <- ",input$cancerSurv, 
                           "genes_list")))
    selectInput("choosen_geneSurv", label=NULL, 
                choices=tmp, selected=tmp[[1]])
  })
  
  #Creating stratum value slider for survival curve plot
  output$plot_help_surv1 <- renderUI({
    if(is.null(input$cancerSurv) || is.null(input$choosen_geneSurv)){
      return(sliderInput("stratum_value_choosen", label="Stratum value: ", 
                  min=0, max=1,
                  value=0.5, step=0.01, round=-2, ticks=FALSE))
    }
    eval(parse(text=paste0("gene <- ",input$cancerSurv, 
                           "data[, '",input$choosen_geneSurv,"']")))
    max <- max(gene, na.rm=TRUE)
    min <- min(gene, na.rm=TRUE)
    med <- median(gene, na.rm=TRUE)
    sliderInput("stratum_value_choosen", label="Stratum value: ", 
                min=round(min,2), max=round(max,2),
                value=med, step=0.01, round=-2, ticks=FALSE)
    
  })
  
  #Creating histogram plot
  output$histogram <- renderPlot({
    if(is.null(input$choosen_geneSurv) || is.null(input$cancerSurv)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
      suppressMessages({
        eval(parse(text=paste0("q <- ",
                               input$cancerSurv, "data[,'",input$choosen_geneSurv,
                               "',drop=FALSE]")))
        q<-ggplot(q,aes_string(names(q)[1]))+ geom_histogram(binwidth = 0.1)+
          scale_y_continuous(name="Quantity")+
          scale_x_continuous(name="mRNA expression value")+
          ggtitle(paste0("mRNA histogram for ",input$choosen_geneSurv))+
          theme(plot.title = element_text(lineheight=20, face="bold"))
        return(q)
      })
    
  })
  
  #Creating survival plot
  output$survplot <- renderPlot({
    if(is.null(input$cancerSurv) || is.null(input$choosen_geneSurv)){
      p <-plot.new()
      text(0.5,0.5,"LOADING...")
      return (p)
    }
      try({
        eval(parse(text=paste0("tmp <- which(colnames(",input$cancerSurv,
                               "data)=='", input$choosen_geneSurv,"')")))
        eval(parse(text=paste0("dane <- ",input$cancerSurv, "data[,c(1,",
                               tmp,",(ncol(",
                               input$cancerSurv,"data)-1):ncol(",
                               input$cancerSurv,"data))]")))
        dane$time <- dane$time/365
        med <- input$stratum_value_choosen
        eval(parse(text=paste0('q<-survfit(Surv(time, status 
                               == "dead")~(', input$choosen_geneSurv, 
                               ">med), data=dane)")))
        suppressMessages({
        q<-autoplot(q,xLab="years",yLab="Survival",
                    title=paste0("Survival plot for mRNA ",
                                 input$choosen_geneSurv,
                                 "\n1-P(death)"),
                    legTitle = "Stratum",
                    legLabs=c(paste0("Patients with value < ",
                                     round(med,4)),
                              paste0("Patients with value > ",
                                     round(med,4))))$plot+scale_y_continuous(
                                       name="", limits=c(0,1)
                                     )+theme(legend.justification=c(0,0),
                                             legend.position=c(0,0))
        return(q)
      })}, silent=TRUE)
  })
  
  #Creating text information under surv - plot
  output$plot_help_surv2 <- renderText({
    if(is.null(input$choosen_geneSurv) || is.null(input$cancerSurv)){
      return ("LOADING...")
    }
    eval(parse(text=paste0("gene <- ",input$cancerSurv, 
                           "data[, '",input$choosen_geneSurv,"']")))
    n_low <- length(gene[gene<=input$stratum_value_choosen])
    n_high <- length(gene[gene>input$stratum_value_choosen])
    
    paste0("Number of observations <=", input$stratum_value_choosen,
           ": ", n_low, br() , 
           "Number of observations >", input$stratum_value_choosen,
           ": ", n_high, br(),
           "Tota number of observations: ", n_low+n_high, br(),br(),
           strong("Attention:"), br(),
           "If the number of observations is small in one of above groups
           the mRNA can seem to be significant, 
           but it does not have to be.",br(),br())
  })
  

   observeEvent(input$save, {
     try({
        suppressMessages({
           eval(parse(text=paste0("tmp <- which(colnames(",input$cancerSurv,
                                  "data)=='", input$choosen_geneSurv,"')")))
           eval(parse(text=paste0("dane <- ",input$cancerSurv, "data[,c(1,",
                                  tmp,",(ncol(",
                                  input$cancerSurv,"data)-1):ncol(",
                                  input$cancerSurv,"data))]")))
           dane$time <- dane$time/365
           med <- input$stratum_value_choosen
           eval(parse(text=paste0('q<-survfit(Surv(time, status 
                                  == "dead")~(', input$choosen_geneSurv, 
                                  ">med), data=dane)")))
           q<-autoplot(q,xLab="years",yLab="Survival",
                       title=paste0("Survival plot for mRNA ",
                                    input$choosen_geneSurv,
                                    "\n1-P(death)"),
                       legTitle = "Stratum",
                       legLabs=c(paste0("Patients with value < ",
                                        round(med,4)),
                                 paste0("Patients with value > ",
                                        round(med,4))))$plot+scale_y_continuous(
                                           name="", limits=c(0,1)
                                        )+theme(legend.justification=c(0,0),
                                                legend.position=c(0,0))
           
           
           lf <- list.files("./www/", pattern=glob2rx("rplot*.jpg"))
           if(length(lf)==0)
           {
              jpeg(sprintf('./www/rplot%s.jpg',1))
              print(q)
              dev.off()
           }
           else
           {
              numbers_of_plots<-unique(na.omit(as.numeric(unlist(strsplit(lf, "[^0-9]+")))))
              directories<-paste0("./www/",lf)
              times<-file.info(directories)$mtime
              which_max_file<-which(times==max(times))
              which_max<-numbers_of_plots[which_max_file]
              if(which_max!=1000)
              {
                 jpeg(sprintf('./www/rplot%s.jpg',which_max+1))
                 print(q)
                 dev.off()
              }
              else
              {
                 jpeg(sprintf('./www/rplot%s.jpg',1))
                 print(q)
                 dev.off()
              }
           }
          
            

            
      
        
        })
     }, silent=TRUE)  
  })

adress<-eventReactive(input$save, {
   
         
         
         lf <- list.files("./www/", pattern=glob2rx("rplot*.jpg"))
         if(length(lf)==0)
         {
            return(list("(rplot1.jpg)","http://mi2.mini.pw.edu.pl:8080/RTCGA/KRWK/shiny/rplot1.jpg"))
         }
         else
         {
            numbers_of_plots<-unique(na.omit(as.numeric(unlist(strsplit(lf, "[^0-9]+")))))
            directories<-paste0("./www/",lf)
            times<-file.info(directories)$mtime
            which_max_file<-which(times==max(times))
            which_max<-numbers_of_plots[which_max_file]
            if(which_max!=1)
            {
               return(list(sprintf("(rplot%s.jpg)",which_max),
                           sprintf("http://mi2.mini.pw.edu.pl:8080/RTCGA/KRWK/shiny/rplot%s.jpg",which_max)))

            }
            else
            {
               return(list("(rplot1.jpg)","http://mi2.mini.pw.edu.pl:8080/RTCGA/KRWK/shiny/rplot1.jpg"))
            }
         }
         
         
         
         
         
         
      })
      output$adres <- renderUI({
         p<-list(adress()[[2]])
         tags$a(paste("Download Kaplan-Meier Plot",adress()[[1]]),href=p,target="_blank")
      })
})

