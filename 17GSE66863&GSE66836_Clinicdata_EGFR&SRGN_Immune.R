library(GSVA)
library(tidyverse)
library(ggpubr)
## Combine driver gene status-----
rm(list=ls())
load("./01Importexpression.Rda")
target<-data.table::fread("./Cluster.txt")
colnames(target)[7]<-"Sample"
table(target$Stage)
target$Stage[str_detect(target$Stage,"III")]<-"III"
target$Stage[target$Stage%in%c("IIa","IIb")]<-"II"
target$Stage[target$Stage%in%c("Ia","Ib")]<-"I"
table(target$Stage)
target$Stage[target$Stage%in%c("I","II")]<-"StageI_II"
target$Stage[target$Stage%in%c("III","IV")]<-"StageIII_IV"
target$age_group[target$Age<65]<-"<65"
target$age_group[target$Age>=65]<-">=65"

### Mutation-----
library(tidyverse)
Targets<-data.table::fread("./Cluster.txt")
egfr<-Targets[,c("Cluster","EGFR_status")]%>%
  arrange(Cluster)%>%
  filter(!is.na(EGFR_status))%>%
  mutate(c_Egfr=paste0(EGFR_status,Cluster))%>%
  add_count(c_Egfr,name="Count")%>%
  distinct(c_Egfr,.keep_all = T)
egfr$EGFR_status<-factor(egfr$EGFR_status,labels  = c("EGFR_mut","EGFR_wt"))
egfr_status<-egfr[,-3]%>%pivot_wider(names_from = EGFR_status,values_from =Count )%>%
  column_to_rownames("Cluster")

kras<-Targets[,c("Cluster","KRAS_status")]%>%
  arrange(Cluster)%>%
  filter(!is.na(KRAS_status))%>%
  mutate(c_Kras=paste0(KRAS_status,Cluster))%>%
  add_count(c_Kras,name="Count")%>%
  distinct(c_Kras,.keep_all = T)

kras_status<-kras[,-3]%>%pivot_wider(names_from = KRAS_status,values_from =Count )%>%
  column_to_rownames("Cluster")
kras$KRAS_status<-factor(kras$KRAS_status,labels  = c("KRAS_mut","KRAS_wt"))
kras_status<-as.table(as.matrix(kras_status))

tp53<-Targets[,c("Cluster","Tp53_status")]%>%
  arrange(Cluster)%>%
  filter(!is.na(Tp53_status))%>%
  mutate(c_Tp53=paste0(Tp53_status,Cluster))%>%
  add_count(c_Tp53,name="Count")%>%
  distinct(c_Tp53,.keep_all = T)
tp53_status<-tp53[,-3]%>%pivot_wider(names_from = Tp53_status,values_from =Count )%>%
  column_to_rownames("Cluster")
tp53$Tp53_status<-factor(tp53$Tp53_status,labels  = c("Tp53_mut","Tp53_wt"))
tp53_status<-as.table(as.matrix(tp53_status))




### Export Clinic----
rm(list=ls())
library(compareGroups)
Clinic<-Targets%>%
  column_to_rownames("title")%>%
  dplyr::select(c("Age","Gender","Stage","Smoking_history","Tp53_status","EGFR_status","KRAS_status","Cluster"))
table(Clinic$Stage)
Clinic$Stage[str_detect(Clinic$Stage,"III")]<-"III"
Clinic$Stage[Clinic$Stage%in%c("IIa","IIb")]<-"II"
Clinic$Stage[Clinic$Stage%in%c("Ia","Ib")]<-"I"
table(Clinic$Age<65)
Clinic$Age<-cut(Clinic$Age,
                breaks=c(min(Clinic$Age,na.rm = T),65,max(Clinic$Age,na.rm = T)),
                labels=c("<65",">=65"),right=F)
table(Clinic$Age)
# <65 >=65 
# 51   60 
restab<-descrTable(Cluster~.,data=Clinic,show.all = T)
export2word(restab,file="./00Clinic_feature.docx")



# Figure S 2A immune-----
rm(list=ls())
target<-data.table::fread("./Cluster.txt")
load("./02Methylation_byChAMP.Rda")
beta0<-beta_norm[,target$geo_accession_m]
rm(beta_norm,Myload,pD)
load("./01Importexpression.Rda")
exp0<-expr_norm%>%column_to_rownames("Gene")%>%as.data.frame() %>% select(target$geo_accession_e)
rm(expr_norm,GSE66863_raw,norm,pda)
dat1<-exp0[c("LCK","CD247","PSTPIP1"),target$geo_accession_e] %>% t() %>% as.data.frame()
dat2<-beta0[c('cg09032544', 'cg07786657', 'cg11683242', 'cg26227523'),target$geo_accession_m] %>% t() %>% as.data.frame()
dat<-cbind(dat1,dat2)
dat<-dat[,c("cg09032544","cg07786657","CD247","cg11683242","LCK","cg26227523","PSTPIP1")]
rownames(dat)<-target$title

load("./02Immunecells_gse66863.rda")
IM<-IM[target$geo_accession_e,]
colnames(IM)
Imm<-IM[,c(2,3,4,19,28,30,21,22,15,16)]
CYT<-exp0[c("GZMA","PRF1"),target$geo_accession_e]%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("geo_accession_e")%>%
  mutate(CYT=1/2*(log(PRF1)+log(GZMA)))
Imm<-cbind(Imm[target$geo_accession_e,],
           CYT=CYT[,"CYT"],
           PDCD1=t(exp0["PDCD1",target$geo_accession_e]))
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




#
# Figure 5 B----
dat1<-beta_norm["cg02851793",Target$geo_accession_m] 
dat2<-exp["SRGN",Target$geo_accession_e] 
dat<-data.frame(cg02851793=dat1,t(dat2),Sample=Target$title,EGFR=Target$EGFR_status) %>% 
  filter(!is.na(EGFR))
table(is.na(dat$EGFR))
library(ggExtra)
library(scales)
library(RColorBrewer)
c("#BEAED4","#E41A1C")
show_col(brewer.pal(8,"Set1"))
color<-brewer.pal(8,"Set1")[1:2]
names(color)<-c("mutated","wild type")
cor.test(dat$cg02851793,dat$SRGN)
pdf("./CorelationofSRGN_methy_revised.pdf",height=4,width=4)
p<-ggplot(dat,aes(x=cg02851793,y=SRGN))+
  geom_point(size = 2,shape = 1, fill = 'gray',aes(color=EGFR))+
  geom_smooth(method=lm,formula= y ~ x)+
  theme_bw()+
  theme(legend.position=c(0.8,0.8))+
  scale_color_manual(values=color)+
  annotate("text", label = paste0("r = -0.52","\n","p = 4.645e-09"),
           x = 0.2, y = 6, size = 4, colour = "#696969",family = "serif", fontface = "italic")+
  theme(axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=12))
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)


dev.off()

with(dat,tapply(SRGN,EGFR,shapiro.test))
with(dat,tapply(cg02851793,EGFR,shapiro.test))
