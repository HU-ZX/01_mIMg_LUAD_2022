rm(list=ls())
library(tidyverse)
load("./Cluster_Phenotype_Plot/gse60644Immunescore.rda")
IM_score<-IM_score%>%column_to_rownames("Sample")
load("01Importdata.rda")
load("02Methylation.Rda")
gse56044_beta<-myNorm
gse56044_exp<-norm[,colnames(myNorm)]
IM_score<-IM_score[colnames(myNorm),]
rm(myNorm,non_norm,norm,pda_co,xxx)
load("/Volumes/T7/Methylation/00Gene_exp_methy_match.rda")
probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe,
                probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
gene<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_u_cis$gene,
               probe_exp_methy_d_trans$gene,probe_exp_methy_u_trans$gene))
probe_exp_methy<-rbind(probe_exp_methy_d_cis,
                       probe_exp_methy_d_trans,
                       probe_exp_methy_u_cis,
                       probe_exp_methy_u_trans)
res_imps<-data.frame()
j=1
for(i in 1:length(probe)){
  b<-gse56044_beta[probe[i],]
  test<-cor.test(IM_score$Im_score,b,method="spearman")
  if(!anyNA(c(test$p.value,test$estimate))){
    if(test$p.value<0.05 & abs(test$estimate)>0.5){
      res_imps[j,1] <- probe[i]
      res_imps[j,2] <- test$estimate
      j=j+1
    } 
  }
}
colnames(res_imps)<-c("Probe","Imps_Correlation")
gse56044_imps<-res_imps
save(gse56044_imps,file="./06Imps_gse56044_sp.rda")

res_imps<-data.frame()
j=1
IM_score<-IM_score[colnames(gse56044_exp),]
for(i in 1:length(gene)){
  b<-as.numeric(gse56044_exp[gene[i],])
  test<-cor.test(IM_score$Im_score,b,method="spearman")
  if(!anyNA(c(test$p.value,test$estimate))){
    if(test$p.value<0.05 & abs(test$estimate)>0.5){
      res_imps[j,1] <- gene[i]
      res_imps[j,2] <- test$estimate
      j=j+1
    } 
  }
}
colnames(res_imps)<-c("Gene","Imps_Correlation")
gse56044_exp_imps<-res_imps
save(gse56044_exp_imps,file="./06Imps_exp_gse56044_sp.rda")
