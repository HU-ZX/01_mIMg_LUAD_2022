library(tidyverse)
library(survival)
library(survminer)
rm(list=ls())
load("/Users/huzixin/Data_Project/Methylation/00mIMg.rda")
load("./02.TPM&Targets.Rda")
target<-data.table::fread("./Cluster_mut.txt")
mIMg<-result
## Multivariate Cox analysis(stage age gender)----
##Gene Multicox----
exp0<-log2(TPM_sym+1)
gene<-unique(mIMg$gene)
rownames(exp0)<-str_replace_all(rownames(exp0),"-","_")
gene<-str_replace_all(gene,"-","_")
Clinic<-target%>%
  filter(!is.na(OS.Time)&!(OS.Time==0))
colnames(Clinic)[c(9,11)]<-c("time","status")
data<-cbind(t(exp0[gene,Clinic$Sample]),
            Clinic[,c(9,11,10,22,23)])

colnames(data)
Surv_ob<-Surv(data$time, data$status)
result<-data.frame(matrix(NA,ncol=5,nrow=length(gene)))
for(i in 1:length(gene)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(c(gene[i],colnames(data)[22:24]),collapse="+"))) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2][1],2)
  PValue <- round(GSum$coefficients[,5][1],3)
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result[i,]<-c(gene[i],HR,names(HR),CI,PValue)}
colnames(result)<-c("gene","HR","group","95%CI","P.val")
result$HR<-as.numeric(result$HR)
result$P.val<-as.numeric(result$P.val)
result1_tcga<-result

save(result1_tcga,file="./tcga_Survival_Results/tcga_multicox1.rda")

##Gene Unicox------
data<-cbind(t(exp0[gene,Clinic$Sample]),
            Clinic[,c(9,11,10)])
Unicox_result<-data.frame(matrix(NA,ncol=4,nrow=length(gene)))
for(i in 1:length(gene)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', gene[i],collapse="+")) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") 
  Unicox_result[i,]<-c(gene[i],HR,CI,PValue)}
colnames(Unicox_result)<-c("gene","HR","95%CI","P.val")
Unicox_result$HR<-as.numeric(Unicox_result$HR) 
Unicox_result$P.val<-as.numeric(Unicox_result$P.val) 
Unicox_result1_tcga<-Unicox_result
save(Unicox_result1_tcga,file="./tcga_Survival_Results/tcga_unicox1.rda")

##Probe Multicox------
load("./04CHAMPline_import.Rda")
beta0<-beta_norm
rm(beta_norm,Myload)
probe<-unique(mIMg$Probe)
Clinic<-target%>%
  filter(!is.na(OS.Time)&!(OS.Time==0))
colnames(Clinic)[c(9,11)]<-c("time","status")

data<-cbind(t(beta0[probe,Clinic$Sample]),
            Clinic[,c(9,11,10,22,23)])

colnames(data)
Surv_ob<-Surv(data$time, data$status)
result<-data.frame(matrix(NA,ncol=5,nrow=length(probe)))
for(i in 1:length(probe)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = probe[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(c(probe[i],colnames(data)[27:29]),collapse="+"))) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2][1],2)
  PValue <- round(GSum$coefficients[,5][1],3)
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result[i,]<-c(probe[i],HR,names(HR),CI,PValue)}
colnames(result)<-c("Probe","HR","group","95%CI","P.val")
result$HR<-as.numeric(result$HR)
result$P.val<-as.numeric(result$P.val)
result1_tcga_probe<-result %>% 
  inner_join(mIMg[,c("gene","Probe")],by="Probe")
save(result1_tcga_probe,file="./tcga_Survival_Results/tcga_multicox1_probe.rda")

## Probe Unicox-----
data<-cbind(t(beta0[probe,Clinic$Sample]),
            Clinic[,c(9,11)])
unicox_result<-data.frame(matrix(NA,ncol=4,nrow=length(probe)))
for(i in 1:length(probe)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = probe[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', probe[i],collapse="+"))
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") 
  unicox_result[i,]<-c(probe[i],HR,CI,PValue)}

colnames(unicox_result)<-c("Probe","HR","95%CI","P.val")
unicox_result$HR<-as.numeric(unicox_result$HR)
unicox_result$P.val<-as.numeric(unicox_result$P.val)
colnames(unicox_result)[2:4]<-paste0(colnames(unicox_result)[2:4],"_unicox")

unicox_result1_tcga_probe<-unicox_result %>% 
  inner_join(mIMg[,c("gene","Probe")],by="Probe")
write.table(unicox_result1_tcga_probe,"./tcga_Survival_Results/02coxprobe.txt",quote=F,sep="\t",col.names=T,row.names=F)                                



