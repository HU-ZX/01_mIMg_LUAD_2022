rm(list=ls())
library(tidyverse)
load("./Cluster_Phenotype_Plot/gse66863Immunescore.rda")
load("./02Methylation_byChAMP.Rda")
load("./01Importexpression.Rda")
rownames(pD)<-pD$geo_accession_m
pD<-pD[colnames(beta_norm),]
identical(colnames(beta_norm),rownames(pD))
colnames(beta_norm)<-pD$title
rownames(pD)<-pD$geo_accession_e
expr_norm<-expr_norm%>%
  column_to_rownames("Gene") %>% as.matrix()
expr_norm<-expr_norm[,pD$geo_accession_e]
colnames(expr_norm)<-pD$title
identical(colnames(expr_norm),colnames(beta_norm))
gse66836_beta<-beta_norm
gse66836_exp<-expr_norm

IM_score<-IM_score %>% column_to_rownames("Sample")
IM_score<-IM_score[pD$geo_accession_e,]
identical(rownames(IM_score),pD$geo_accession_e)
rownames(IM_score)<-pD$title


load("/Volumes/T7/Methylation/00Gene_exp_methy_match.rda")
probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe,
                probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
gene<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_u_cis$gene,
               probe_exp_methy_d_trans$gene,probe_exp_methy_u_trans$gene))
probe_exp_methy<-rbind(probe_exp_methy_d_cis,
                       probe_exp_methy_d_trans,
                       probe_exp_methy_u_cis,
                       probe_exp_methy_u_trans)
rm(beta_norm,Myload,GSE66863_raw,expr_norm,pda,norm)
## Probe-----
IM_score<-IM_score[colnames(gse66836_beta),]
res_imps<-data.frame()
j=1
for(i in 1:length(probe)){
  b<-gse66836_beta[probe[i],]
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
gse66836_imps<-res_imps
save(gse66836_imps,file="./06Imps_gse66836_sp.rda")

## Gene----
res_imps<-data.frame()
j=1
IM_score<-IM_score[colnames(gse66836_exp),]
for(i in 1:length(gene)){
  b<-as.numeric(gse66836_exp[gene[i],])
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
gse66836_exp_imps<-res_imps
save(gse66836_exp_imps,file="./06Imps_exp_gse66836_sp.rda")


