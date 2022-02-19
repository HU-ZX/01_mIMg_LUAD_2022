rm(list=ls())
library(tidyverse)
load("02Cluster_IM.Rda")
load("./04CHAMPline_import.Rda")
load("/Volumes/T7/Methylation/00Gene_exp_methy_match.rda")
probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe,
                probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
gene<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_u_cis$gene,
               probe_exp_methy_d_trans$gene,probe_exp_methy_u_trans$gene))
probe_exp_methy<-rbind(probe_exp_methy_d_cis,
                       probe_exp_methy_d_trans,
                       probe_exp_methy_u_cis,
                       probe_exp_methy_u_trans)
load("02.TPM&Targets.Rda")
tcga_exp<-log2(TPM_sym+1)%>%as.data.frame()
Targets<-data.table::fread("./Cluster.txt")%>%arrange(Cluster)
tcga_beta<-beta_norm
IM_score<-IM%>%mutate(Im_score=Effector+Suppressor+MHC+Th1+Th2+Immunoinhibitor+Immunostimulator)%>%
  rownames_to_column("Sample")%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")%>%
  column_to_rownames("Sample")
CYT<-CYT%>%column_to_rownames("Sample")
rm(beta_norm,Myload,res_gsva,Type,IM,IPS,TPM,TPM_sym)

IM_score<-IM_score[colnames(tcga_beta),]
res_imps<-data.frame()
j=1
for(i in 1:length(probe)){
  b<-tcga_beta[probe[i],]
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
tcga_imps<-res_imps
save(tcga_imps,file="./06Imps_tcga_sp.rda")

res_imps<-data.frame()
j=1
IM_score<-IM_score[colnames(tcga_exp),]
for(i in 1:length(gene)){
  b<-as.numeric(tcga_exp[gene[i],])
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
tcga_exp_imps<-res_imps
save(tcga_exp_imps,file="./06Imps_exp_tcga_sp.rda")
cor.test(as.numeric(IM_score[Targets$Sample,"Im_score"]), as.numeric(CYT[Targets$Sample,"CYT"]))

