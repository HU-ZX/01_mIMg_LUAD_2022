###Export Clinic-----
library(tidyverse)
Targets<-data.table::fread("./Cluster.txt")
library(compareGroups)
Clinic<-Targets%>%
  column_to_rownames("title")%>%
  dplyr::select(c("age","sex","Cluster"))%>%
  rename("Gender"="sex",
         "Age"="age")

table(Clinic$Age<65)
Clinic$Age<-cut(Clinic$Age,
                breaks=c(min(Clinic$Age,na.rm = T),65,max(Clinic$Age,na.rm = T)),
                labels=c("<65",">=65"),right=F)

table(Clinic$Age)
restab<-descrTable(Cluster~.,data=Clinic,show.all = T)
export2csv(restab,file="./00Clinic_feature.csv")

### Figure S2B Immune------
rm(list=ls())
load("01Importdata.rda")
target<-data.table::fread("./Cluster.txt")%>%arrange(Cluster)
exp0<-norm[,target$Sample]
rm(non_norm,norm,pda_co,xxx)
load("02Methylation.Rda")
beta0<-myNorm[,target$Sample]
rm(myNorm,pD)
dat1<-exp0[c("LCK","CD247","PSTPIP1"),target$Sample] %>% t() %>% as.data.frame()
dat2<-beta0[c('cg09032544', 'cg07786657', 'cg11683242', 'cg26227523'),target$Sample] %>% t() %>% as.data.frame()
dat<-cbind(dat1,dat2)
dat<-dat[,c("cg09032544","cg07786657","CD247","cg11683242","LCK","cg26227523","PSTPIP1")]
load("./02Cluster.Rda")
IM<-t(res_gsva)[target$Sample,c(2,3,4,19,28,30,21,22,15,16)]
rm(res_gsva)
CYT<-exp0[c("GZMA","PRF1"),]%>%
  t()%>%as.data.frame()%>%
  mutate(CYT=1/2*(log(PRF1)+log(GZMA)))
Imm<-cbind(IM[target$Sample,],
           CYT=CYT[target$Sample,"CYT"],
           PDCD1=t(exp0["PDCD1",target$Sample]))
IM<-IM[target$Sample,]
colnames(IM)
rm(CYT,IM)
library(corrplot)
corr<-WGCNA::cor(dat,Imm, use = "p",method="spearman")
nSamples<-nrow(dat)
p.mat<-WGCNA::corPvalueStudent(corr, nSamples)
pdf("./methy_Im_Corr.pdf",width=10,height=10)
corrplot(corr,p.mat=p.mat,sig.level=0.05,addshade = "all",
         addCoef.col = 'grey',col=pal_gsea("default", n = 100, reverse = F)(100),
         method = "circle",number.cex=1.2,tl.col="black",tl.cex=1)
dev.off()



