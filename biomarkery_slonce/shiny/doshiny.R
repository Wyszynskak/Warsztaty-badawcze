# library(ggplot2)
# library(survMisc)
# library(RTCGA.mRNA)
# library(RTCGA.clinical)
# cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG",
#                 "LUAD","LUSC","OV","READ","UCEC")
# setwd("C:/Users/wyszy/OneDrive/Studia/Magisterka/Rok 2 sem 1/Warsztaty badawcze/Projekt końcowy/biomarkery_slonce")
# for (j in 1:length(cancer_types))
# {
#    p<-sprintf("%s.mRNA",cancer_types[j])
#    q<-sprintf("%s.clinical",cancer_types[j])
#    namep<-data(list=p)
#    nameq<-data(list=q)
#    
#    
#    mBRCA<-.GlobalEnv[[namep]]
#    cBRCA<-.GlobalEnv[[nameq]]
#    cBRCA <- data.frame(time1=as.numeric(as.character(cBRCA$patient.days_to_death)),
#                        time2=as.numeric(as.character(cBRCA$patient.days_to_last_followup)),
#                        status = cBRCA$patient.vital_status,
#                        barcode = cBRCA$patient.bcr_patient_barcode,
#                        race=cBRCA$patient.race,
#                        therapy=cBRCA$patient.drugs.drug.therapy_types.therapy_type)
#    cBRCA$time <- ifelse(is.na(cBRCA$time1), cBRCA$time2, cBRCA$time1)
#    mBRCA$bcr_patient_barcode<-tolower(substr(mBRCA$bcr_patient_barcode, 1, 12))
#    allBRCA<-merge(mBRCA,cBRCA,by.x="bcr_patient_barcode",by.y="barcode")
#    length_gen<-length(2:ncol(mBRCA))
#    histy<-list()
#    survy<-list()
#    geny<-c()
#    pvalues <- c()
#    for(i in 2:ncol(mBRCA))
#    {
#       tylko_jeden_gen <- allBRCA[,c(1,i,(ncol(mBRCA)+1):ncol(allBRCA))]
#       geny[i-1]<- names(tylko_jeden_gen)[2]
#       tylko_jeden_gen <- tylko_jeden_gen[,c(1,2,5,8)]
#       names(tylko_jeden_gen)[2]<-"val"
#       med<-median(na.omit(tylko_jeden_gen$val))
#       prog <- 0.2*length(tylko_jeden_gen$val)
#       tryCatch({
#          if(length(tylko_jeden_gen$val[tylko_jeden_gen$val>med])>prog & length(
#             tylko_jeden_gen$val[tylko_jeden_gen$val<med])>prog){
#             model <- survdiff(Surv(time, status == "dead")~(val>med),
#                               data=tylko_jeden_gen)
#             pvalues[i-1] <- pchisq(model$chisq, 1, lower.tail=FALSE)
#          } else{
#             pvalues[i-1] <- 100; #Wystarczajaco duza wartosc zeby zmienna byla 
#                                  #"nieistotna"
#          }
#         names(pvalues)[i-1] <- geny[i-1]
#       },
#       error = function(e) {
#          print(paste0("zlo_model",i))
#       })
# 
# #          tryCatch({
# #             p<-survfit(Surv(time, status == "dead")~(val>med), data=tylko_jeden_gen)
# #             survy[[i-1]]<-p
# #             names(survy)[i-1] <- geny[i-1]
# #             
# #          },
# #          error = function(e) {
# #             print(paste0("zlo_surv",i))
# #          })
#          tryCatch({
# 
#             q <- tylko_jeden_gen[,'val',drop=FALSE]
#             typeof(q)
#             
#             histy[[i-1]]<-q
#             names(histy)[i-1] <- geny[i-1]
#          },
#          error = function(e) {
#             print(paste0("zlo_hist",i))
#          })
#       if(i/1000==round(i/1000))
#       {
#          print(paste(j,i))
#       }   
#    }
# 
# #    #### Wybranie 1000 najbardziej istotnych genow
# #    names(pvalues) <- geny
# #    pvalues <- sort(pvalues)[1:1000]
# #    geny_interesujace <- names(pvalues)
# #    #save(geny_interesujace, file="D:/geny_interesujace_BRCA.RData")
# #    survy <- survy[geny_interesujace]
# #    histy <- histy[geny_interesujace]
#    
# #    if(!(dir.exists(paste0("./",cancer_types[j]))))
# #    {
# #       dir.create(paste0("./",cancer_types[j])) 
# #    }
#    #saveRDS(geny_interesujace,paste0("./shiny/data/",cancer_types[j],"geny.RDS"))
#    #saveRDS(survy,paste0("./",cancer_types[j],"/survy.RDS"))
#    pvalues <- sort(pvalues)
#    saveRDS(histy,paste0("./shiny/data/",cancer_types[j],"dane.RDS"))
#    saveRDS(pvalues,paste0("./shiny/data/",cancer_types[j],"pvalue.RDS"))
# }



##############################################################################
###############################################################################
###############################################################################
#Podejcie drugie
library(ggplot2)
library(survMisc)
library(RTCGA.mRNA)
library(RTCGA.clinical)
cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG",
                "LUAD","LUSC","OV","READ","UCEC")
setwd("C:/Users/wyszy/OneDrive/Studia/Magisterka/Rok 2 sem 1/Warsztaty badawcze/Projekt końcowy/biomarkery_slonce")

for (j in 1:length(cancer_types))
{
  p<-sprintf("%s.mRNA",cancer_types[j])
  q<-sprintf("%s.clinical",cancer_types[j])
  namep<-data(list=p)
  nameq<-data(list=q)
  
  
  mBRCA<-.GlobalEnv[[namep]]
  cBRCA<-.GlobalEnv[[nameq]]
  cBRCA <- data.frame(time1=as.numeric(as.character(cBRCA$patient.days_to_death)),
                      time2=as.numeric(as.character(cBRCA$patient.days_to_last_followup)),
                      status = cBRCA$patient.vital_status,
                      barcode = cBRCA$patient.bcr_patient_barcode,
                      race=cBRCA$patient.race,
                      therapy=cBRCA$patient.drugs.drug.therapy_types.therapy_type)
  cBRCA$time <- ifelse(is.na(cBRCA$time1), cBRCA$time2, cBRCA$time1)
  mBRCA$bcr_patient_barcode<-tolower(substr(mBRCA$bcr_patient_barcode, 1, 12))
  cBRCA <- cBRCA[ , c(3,4,7)]
  allBRCA<-merge(mBRCA,cBRCA,by.x="bcr_patient_barcode",by.y="barcode")
  length_gen<-length(2:ncol(mBRCA))
  pvalues <- c()
  pvalues_num <- c()
  for(i in 2:ncol(mBRCA))
  {
    tylko_jeden_gen <- allBRCA[,c(1,i,(ncol(mBRCA)+1):ncol(allBRCA))]
    tmp <- names(tylko_jeden_gen)[2]
    names(tylko_jeden_gen)[2]<-"val"
    med<-median(na.omit(tylko_jeden_gen$val))
    prog <- 0.2*length(tylko_jeden_gen$val)
    tryCatch({
      if(length(tylko_jeden_gen$val[tylko_jeden_gen$val>med])>prog & length(
        tylko_jeden_gen$val[tylko_jeden_gen$val<med])>prog){
        model <- survdiff(Surv(time, status == "dead")~(val>med),
                          data=tylko_jeden_gen)
        pvalues_num[i-1] <- pchisq(model$chisq, 1, lower.tail=FALSE)
        pvalues[i-1]<-as.character(signif(pvalues_num[i-1],1))
      } else{
        pvalues[i-1] <- "100"; #Wystarczajaco duza wartosc zeby zmienna byla 
        #"nieistotna"
      }
      names(pvalues)[i-1] <- tmp
    },
    error = function(e) {
      print(paste0("zlo_model",i))
    })
    if(i/1000==round(i/1000))
    {
      print(paste(j,i))
    }   
  }
  pvalues <- sort(pvalues)
  text1 <- paste0('"',names(pvalues), "\n (p-value: ", 
                  ifelse(pvalues<100, format(pvalues,scientific=FALSE),
                         "-"),')"="',names(pvalues),'"', collapse=", ")
  eval(parse(text=paste0("genes_list <- list(", text1, ")")))
  saveRDS(allBRCA,paste0("./shiny/data/",cancer_types[j],"dane.RDS"))
  saveRDS(genes_list,paste0("./shiny/data/",cancer_types[j],"genes_list.RDS"))
}


