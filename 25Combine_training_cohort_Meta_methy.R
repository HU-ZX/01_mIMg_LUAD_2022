library(tidyverse)
library(ggpubr)
rm(list=ls())
### Cluster3 vs Cluster2
## Improt myDMP-----
TCGA_Methy<-data.table::fread("./Download/TCGA_MyDMP.txt")
GSE56044_Methy<-data.table::fread("./GEO/GSE60644_GSE60645_Data/GSE56044_MyDMP.txt")
GSE66836_Methy<-data.table::fread("./GEO/GSE66863_GSE66836_Data/GSE66836_myDMP.txt")
methy_ls<-list(TCGA_Methy,GSE56044_Methy,GSE66836_Methy)
## 01 intersected DMP----
## Probe-----

Up_probe<-Reduce(intersect,list(TCGA_Methy$Probe[TCGA_Methy$logFC > 0],
                                GSE56044_Methy$Probe[GSE56044_Methy$logFC > 0],
                                GSE66836_Methy$Probe[GSE66836_Methy$logFC > 0]))

Up<-lapply(methy_ls,function(x){
  x%>%filter(Probe%in%Up_probe)%>%return()
})

Down_probe<-Reduce(intersect,list(TCGA_Methy$Probe[TCGA_Methy$logFC < 0],
                                  GSE56044_Methy$Probe[GSE56044_Methy$logFC < 0],
                                  GSE66836_Methy$Probe[GSE66836_Methy$logFC < 0]))
Down<-lapply(methy_ls,function(x){
  x%>%filter(Probe%in%Down_probe)%>%return()
})

probe<-unique(c(Up_probe,Down_probe))
save(probe,file = "./00Probe_match.rda")
## 02 co_expression methylation------
rm(list = ls())
setwd("/Users/huzixin/Data_Project/Methylation/")
load("./00Probe_match.rda")
load("./GEO/GSE60644_GSE60645_Data/05Diffgene_gse60644.rda")
load("./GEO/GSE66863_GSE66836_Data/05Diffgene_gse66863.rda")
load("./Download/05Diffgene_tcga.rda")
load("./GEO/GSE60644_GSE60645_Data/05exp_methy_logFC.rda")
load("./GEO/GSE66863_GSE66836_Data/05exp_methy_logFC.rda")
load("./Download/05exp_methy_logFC.rda")

probe_exp_methy_u_cis<-Reduce(intersect,list(tcga_Expr_Methy$Probe[tcga_Expr_Methy$logFC_m>0&tcga_Expr_Methy$logFC_g>0],
                                             gse60644_Expr_Methy$Probe[gse60644_Expr_Methy$logFC_m>0&gse60644_Expr_Methy$logFC_g>0],
                                             gse66836_Expr_Methy$Probe[gse66836_Expr_Methy$logFC_m>0&gse66836_Expr_Methy$logFC_g>0]))


probe_exp_methy_d_cis<-Reduce(intersect,list(tcga_Expr_Methy$Probe[tcga_Expr_Methy$logFC_m<0&tcga_Expr_Methy$logFC_g<0],
                                             gse60644_Expr_Methy$Probe[gse60644_Expr_Methy$logFC_m<0&gse60644_Expr_Methy$logFC_g<0],
                                             gse66836_Expr_Methy$Probe[gse66836_Expr_Methy$logFC_m<0&gse66836_Expr_Methy$logFC_g<0]))

probe_exp_methy_u_trans<-Reduce(intersect,list(tcga_Expr_Methy$Probe[tcga_Expr_Methy$logFC_m>0&tcga_Expr_Methy$logFC_g<0],
                                               gse60644_Expr_Methy$Probe[gse60644_Expr_Methy$logFC_m>0&gse60644_Expr_Methy$logFC_g<0],
                                               gse66836_Expr_Methy$Probe[gse66836_Expr_Methy$logFC_m>0&gse66836_Expr_Methy$logFC_g<0]))

probe_exp_methy_d_trans<-Reduce(intersect,list(tcga_Expr_Methy$Probe[tcga_Expr_Methy$logFC_m<0&tcga_Expr_Methy$logFC_g>0],
                                               gse60644_Expr_Methy$Probe[gse60644_Expr_Methy$logFC_m<0&gse60644_Expr_Methy$logFC_g>0],
                                               gse66836_Expr_Methy$Probe[gse66836_Expr_Methy$logFC_m<0&gse66836_Expr_Methy$logFC_g>0]))


colnames(probe_exp_methy_d_cis)
colnames(tcga_Expr_Methy)[3:4]<-paste0(colnames(tcga_Expr_Methy)[3:4],"_tcga")
colnames(gse60644_Expr_Methy)[3:4]<-paste0(colnames(gse60644_Expr_Methy)[3:4],"_gse60644")
colnames(gse66836_Expr_Methy)[3:4]<-paste0(colnames(gse66836_Expr_Methy)[3:4],"_gse66836")
probe_exp_methy_d_cis<-tcga_Expr_Methy%>%
  filter(Probe%in%probe_exp_methy_d_cis)%>%
  inner_join(gse60644_Expr_Methy[gse60644_Expr_Methy$Probe%in%probe_exp_methy_d_cis,-2],by="Probe")%>%
  inner_join(gse66836_Expr_Methy[gse66836_Expr_Methy$Probe%in%probe_exp_methy_d_cis,-2],by="Probe")%>%
  dplyr::select( "Probe","gene","logFC_m_tcga","logFC_g_tcga",
                 "logFC_m_gse60644","logFC_g_gse60644",
                 "logFC_m_gse66836","logFC_g_gse66836",
                 everything())
probe_exp_methy_d_trans<-tcga_Expr_Methy%>%
  filter(Probe%in%probe_exp_methy_d_trans)%>%
  inner_join(gse60644_Expr_Methy[gse60644_Expr_Methy$Probe%in%probe_exp_methy_d_trans,-2],by="Probe")%>%
  inner_join(gse66836_Expr_Methy[gse66836_Expr_Methy$Probe%in%probe_exp_methy_d_trans,-2],by="Probe")%>%
  dplyr::select( "Probe","gene","logFC_m_tcga","logFC_g_tcga",
                 "logFC_m_gse60644","logFC_g_gse60644",
                 "logFC_m_gse66836","logFC_g_gse66836",
                 everything())
probe_exp_methy_u_trans<-tcga_Expr_Methy%>%
  filter(Probe%in%probe_exp_methy_u_trans)%>%
  inner_join(gse60644_Expr_Methy[gse60644_Expr_Methy$Probe%in%probe_exp_methy_u_trans,-2],by="Probe")%>%
  inner_join(gse66836_Expr_Methy[gse66836_Expr_Methy$Probe%in%probe_exp_methy_u_trans,-2],by="Probe")%>%
  dplyr::select( "Probe","gene","logFC_m_tcga","logFC_g_tcga",
                 "logFC_m_gse60644","logFC_g_gse60644",
                 "logFC_m_gse66836","logFC_g_gse66836",
                 everything())
probe_exp_methy_u_cis<-tcga_Expr_Methy%>%
  filter(Probe%in%probe_exp_methy_u_cis)%>%
  inner_join(gse60644_Expr_Methy[gse60644_Expr_Methy$Probe%in%probe_exp_methy_u_cis,-2],by="Probe")%>%
  inner_join(gse66836_Expr_Methy[gse66836_Expr_Methy$Probe%in%probe_exp_methy_u_cis,-2],by="Probe")%>%
  dplyr::select( "Probe","gene","logFC_m_tcga","logFC_g_tcga",
                 "logFC_m_gse60644","logFC_g_gse60644",
                 "logFC_m_gse66836","logFC_g_gse66836",
                 everything())
save(probe_exp_methy_d_cis,probe_exp_methy_u_trans,probe_exp_methy_d_trans,probe_exp_methy_u_cis,file="./00Gene_exp_methy_match.rda")
## 03 combination correlation between methylation and expression----
rm(list=ls())
library(tidyverse)
load("./00Gene_exp_methy_match.rda")
load("./GEO/GSE60644_GSE60645_Data/04Co_exp_methy_probe_fl_sp.rda")
load("./GEO/GSE66863_GSE66836_Data/04Co_exp_methy_probe_fl_sp.rda")
load("./Download/04Co_exp_methy_probe_fl_sp.rda")
probe_fl<-probe_fl_tcga%>%
  inner_join(probe_fl_56044[,1:3],by="Probe")%>%
  inner_join(probe_fl_66836[,1:3],by="Probe")
probe_exp_meth<-rbind(probe_exp_methy_d_cis,
                      probe_exp_methy_d_trans,
                      probe_exp_methy_u_cis,
                      probe_exp_methy_u_trans
                      
)
probe_exp_meth_cr<-probe_fl%>%
  inner_join(probe_exp_meth,by="Probe")
save(probe_exp_meth_cr,file="./00Probe_filtercr_sp.rda")

load("./00Probe_filtercr_sp.rda")
colnames(probe_exp_meth_cr)
probe_exp_meth_cr<-probe_exp_meth_cr[,c(1,9,4,
                                        16:22,2:3,10:11,7:8,14:15,5:6,12:13)]

write.csv(probe_exp_meth_cr,file="./TableS5probe_exp_meth_cr.csv")






