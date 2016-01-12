library(ggplot2)
library(survMisc)
library(RTCGA.mRNA)
library(RTCGA.clinical)
cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG","LUAD","LUSC","OV","READ","UCEC")
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
   allBRCA<-merge(mBRCA,cBRCA,by.x="bcr_patient_barcode",by.y="barcode")
   nrow(allBRCA)
   statsy<-c()
   for(i in 2:ncol(mBRCA))
   {
      tylko_jeden_gen <- allBRCA[,c(1,i,(ncol(mBRCA)+1):ncol(allBRCA))]
      names(tylko_jeden_gen)[2]<-"val"
      
      med<-median(na.omit(tylko_jeden_gen$val))
      len_glob_val<-length(tylko_jeden_gen$val)
      len_above_val<-sum(na.omit(tylko_jeden_gen$val)>med)
      len_below_val<-sum(na.omit(tylko_jeden_gen$val)<med)
      if(!(len_above_val/len_glob_val<0.2|len_below_val/len_glob_val<0.2))
      {
         tryCatch({
            p<-survdiff(Surv(time,status == "dead")~(val>med), data=tylko_jeden_gen)
            p_val<-1-pchisq(p$chisq,1)
         },
         error = function(e) {
            p_val<-2
            print(paste0("zlo",i))
         })
      }
      else
      {
         p_val<-2
      }
      statsy<-c(statsy,p_val)
   }
   names(statsy)<-names(mBRCA)[2:ncol(mBRCA)]
   if(!(dir.exists(paste0("./",cancer_types[j]))))
   {
      dir.create(paste0("./",cancer_types[j])) 
   }
   saveRDS(statsy,paste0("./",cancer_types[j],"/p-values.RDS"))
   dat_min<-allBRCA[,c(1,which(statsy==min(statsy)),(ncol(mBRCA)+1):ncol(allBRCA))]
   title<-names(dat_min)[2]
   names(dat_min)[2]<-"val"
   med<-median(na.omit(dat_min$val))
   ob1 <- survfit(Surv(time, status == "dead")~(val>med), data=dat_min)
   p<-autoplot(ob1,title = title)
   png(filename=paste0("./",cancer_types[j],"/",title,".png"))
   print(p)
   dev.off()
   rm(mBRCA,cBRCA,tylko_jeden_gen,statsy)
}

