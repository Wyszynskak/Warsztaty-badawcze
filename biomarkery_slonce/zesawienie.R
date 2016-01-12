cancer_types<-c("BRCA","COAD","COADREAD","GBMLGG","KIPAN","KIRC","KIRP","LGG","LUAD","LUSC","OV","READ","UCEC")

most_good_cancer<-list()
for(cancer in cancer_types)
{
   p<-readRDS(paste0("./",cancer,"/p-values.RDS"))
   p<-sort(p,decreasing=)
   p1<-p[1:100]
   if(sum(names(p1)=='SNIP1'))
   print(cancer)
   most_good_cancer[[which(cancer_types==cancer)]]<-p1
}
most_good<-most_good_cancer[[1]]
for(i in 2:length(cancer_types))
{
  most_good<-c(most_good,most_good_cancer[[i]])   
}   
q<-sort(table(names(most_good))) 
statsy<-readRDS(paste0("./",'COAD',"/p-values.RDS"))

dat_min<-allBRCA[,c(1,(which(names(statsy)=='SNIP1')+1),(ncol(mBRCA)+1):ncol(allBRCA))]
title<-names(dat_min)[2]
names(dat_min)[2]<-"val"
med<-median(na.omit(dat_min$val))
ob1 <- survfit(Surv(time, status == "dead")~(val>med), data=dat_min)
autoplot(ob1,title = title)
head(dat_min)
