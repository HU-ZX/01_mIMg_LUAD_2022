# GSE11969----
## Data import----

GSE11969_sm<-getGEO(GEO = "GSE11969", 
                    filename = "./GSE11969_series_matrix.txt",getGPL = F)
GSE11969_pd<-pData(GSE11969_sm)
colnames(GSE11969_pd)
target<-GSE11969_pd[,c(1,2,10,12:25)]
target[1,]

colnames(target)<-c("title","Sample","Annotation",
                    "Age","Sex","Histology","Smoking",
                    "TNM","Stage","status","time",
                    "Evidence_of_Relapse","Site_of_Relapse",
                    "EGFR","Kras","p53","gefitinib")
table(target$Histology)
target<-target%>%
  mutate(Annotation=str_replace(sapply(str_split(Annotation,":"),"[",2)," ",""),
         Age=as.numeric(sapply(str_split(Age,":"),"[",2)),
         Sex=str_replace(sapply(str_split(Sex,":"),"[",2)," ",""),
         Histology=str_replace(sapply(str_split(Histology,":"),"[",2)," ",""),
         Smoking=as.numeric(sapply(str_split(Smoking,":"),"[",2)),
         TNM=str_replace(sapply(str_split(TNM,":"),"[",2)," ",""),
         Stage=str_replace(sapply(str_split(Stage,":"),"[",2)," ",""),
         status=str_replace(sapply(str_split(status,":"),"[",2)," ",""),
         time=as.numeric(sapply(str_split(time,":"),"[",2)),
         Evidence_of_Relapse=str_replace(sapply(str_split(Evidence_of_Relapse,":"),"[",2)," ",""),
         Site_of_Relapse=str_replace(sapply(str_split(Site_of_Relapse,":"),"[",2)," ",""),
         EGFR=str_replace(sapply(str_split(EGFR,":"),"[",2)," ",""),
         Kras=str_replace(sapply(str_split(Kras,":"),"[",2)," ",""),
         p53=str_replace(sapply(str_split(p53,":"),"[",2)," ",""),
         gefitinib=str_replace(sapply(str_split(gefitinib,":"),"[",2)," ",""))%>%
  filter(Histology=="AD")%>%
  mutate(T_stage=substring(TNM,1,2),
         N_stage=substring(TNM,3,4),
         M_stage=substring(TNM,5,6))
target$status<-as.numeric(ifelse(target$status=="Alive","0","1"))
files<-data.frame(FileName=list.files("./GSE11969_RAW"),
                  Sample=str_extract(list.files("./GSE11969_RAW/"),"GSM\\d+"))
target<-target%>%inner_join(files,by="Sample")
rownames(target)<-target$Sample
GSE11969_expr<-exprs(GSE11969_sm)[,target$Sample]
GSE11969_sm<-ExpressionSet(assayData = GSE11969_expr, 
                           phenoData = AnnotatedDataFrame(data = target),
                           annotation = GSE11969_sm@annotation)
save(GSE11969_sm,file="./01GSE11969_Series_matrix.rda")

rm(list=ls())
load("./01GSE11969_Series_matrix.rda")
library(limma)
target<-pData(GSE11969_sm)
GSE11969_raw <- read.maimages(files = target$FileName,
                              source = "agilent",
                              path = "./GSE11969_RAW",
                              names = target$Sample,
                              other.columns = "gIsWellAboveBG",
                              green.only = T)
## probe type summary
table(GSE11969_raw$genes$ControlType)
head(GSE11969_raw$E)[, 1:5]
dim(GSE11969_raw)
## add targets info
GSE11969_raw$targets <- target

eset <- ExpressionSet(assayData = GSE11969_raw$E, 
                      phenoData = AnnotatedDataFrame(data = GSE11969_raw$targets))
library(arrayQualityMetrics)
arrayQualityMetrics(eset, 
                    outdir = "./02GSE11969_QC", 
                    do.logtransform = TRUE,
                    force = TRUE)
dev.off()
bgc <-  backgroundCorrect(RG = GSE11969_raw, 
                          method = "normexp",
                          offset = 50,
                          normexp.method = "mle")

norm <- normalizeBetweenArrays(bgc, method = "quantile")

norm$targets<-target
rownames(norm$E)<-norm$genes$ProbeName
gse11969_norm<-norm

expr_norm <- norm$E

sum(is.na(expr_norm))
boxplot(expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
gpl6480_anno<-read_tsv("./GPL6480-9577.txt", comment = "#")%>%
  select("ID","GENE_SYMBOL","GENE")%>%
  filter(!(str_detect(GENE_SYMBOL,"\\/")))%>%
  rename(Gene_Symbol=GENE_SYMBOL)
vali_exp<-expr_norm %>%as.data.frame()%>%
  rownames_to_column("ID")%>%
  inner_join(gpl6480_anno[,1:2],by="ID")%>%
  dplyr::rename("gene"="Gene_Symbol")%>%
  dplyr::select(-"ID")%>%
  aggregate(.~ gene, data = ., max)%>%
  column_to_rownames("gene")
gse11969_exp<-vali_exp
save(gse11969_exp,gse11969_norm,file="./01GSE11969_90samples.rda")

# GSE13213----
## Data import--------
rm(list=ls())
library(GEOquery)
library(tidyverse)
GSE13213_sm<-getGEO(GEO = "GSE13213", 
                    filename = "./GSE13213_series_matrix.txt",getGPL = F)
GSE13213_pd<-pData(GSE13213_sm)
colnames(GSE13213_pd)
target<-GSE13213_pd[,c(1,2,10:24)]
target[1,]
colnames(target)<-c("title","Sample","Annotation","Cohort",
                    "Age","Sex","Histology","Smoking",
                    "TNM","Stage","status","time",
                    "Evidence_of_Relapse","Site_of_Relapse",
                    "EGFR","Kras","p53")
target<-target%>%
  mutate(Annotation=str_replace(sapply(str_split(Annotation,":"),"[",2)," ",""),
         Cohort=str_replace(sapply(str_split(Cohort,":"),"[",2)," ",""),
         Age=as.numeric(sapply(str_split(Age,":"),"[",2)),
         Sex=str_replace(sapply(str_split(Sex,":"),"[",2)," ",""),
         Histology=str_replace(sapply(str_split(Histology,":"),"[",2)," ",""),
         Smoking=as.numeric(sapply(str_split(Smoking,":"),"[",2)),
         TNM=str_replace(sapply(str_split(TNM,":"),"[",2)," ",""),
         Stage=str_replace(sapply(str_split(Stage,":"),"[",2)," ",""),
         status=str_replace(sapply(str_split(status,":"),"[",2)," ",""),
         time=as.numeric(sapply(str_split(time,":"),"[",2)),
         Evidence_of_Relapse=str_replace(sapply(str_split(Evidence_of_Relapse,":"),"[",2)," ",""),
         Site_of_Relapse=str_replace(sapply(str_split(Site_of_Relapse,":"),"[",2)," ",""),
         EGFR=str_replace(sapply(str_split(EGFR,":"),"[",2)," ",""),
         Kras=str_replace(sapply(str_split(Kras,":"),"[",2)," ",""),
         p53=str_replace(sapply(str_split(p53,":"),"[",2)," ",""))%>%
  mutate(T_stage=substring(TNM,1,2),
         N_stage=substring(TNM,3,4),
         M_stage=substring(TNM,5,6))
target$status<-as.numeric(ifelse(target$status=="Alive","0","1"))

files<-data.frame(FileName=list.files("./GSE13213_RAW/"),
                  Sample=str_extract(list.files("./GSE13213_RAW/"),"GSM\\d+"))
target<-target%>%inner_join(files,by="Sample")
rownames(target)<-target$Sample
GSE13213_expr<-exprs(GSE13213_sm)[,target$Sample]
GSE13213_sm<-ExpressionSet(assayData = GSE13213_expr, 
                           phenoData = AnnotatedDataFrame(data = target),
                           annotation = GSE13213_sm@annotation)
save(GSE13213_sm,file="./01GSE13213_Series_matrix.rda")

rm(list=ls())
load("./01GSE13213_Series_matrix.rda")
library(limma)
target<-pData(GSE13213_sm)
GSE13213_raw <- read.maimages(files = target$FileName,
                              source = "agilent",
                              path = "./GSE13213_RAW",
                              names = target$Sample,
                              other.columns = "gIsWellAboveBG",
                              green.only = T)

## probe type summary
table(GSE13213_raw$genes$ControlType)
head(GSE13213_raw$E)[, 1:5]
dim(GSE13213_raw)
## add targets info
GSE13213_raw$targets <- target

eset <- ExpressionSet(assayData = GSE13213_raw$E, 
                      phenoData = AnnotatedDataFrame(data = GSE13213_raw$targets))
library(arrayQualityMetrics)
arrayQualityMetrics(eset, 
                    outdir = "./02GSE13213_QC_raw", 
                    do.logtransform = TRUE,
                    force = TRUE)
dev.off()



bgc <-  backgroundCorrect(RG = GSE13213_raw, 
                          method = "normexp",
                          offset = 50,
                          normexp.method = "mle")

norm <- normalizeBetweenArrays(bgc, method = "quantile")

norm$targets<-target
rownames(norm$E)<-norm$genes$ProbeName
gse13213_norm<-norm

expr_norm <- norm$E

sum(is.na(expr_norm))
boxplot(expr_norm,
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
gpl6480_anno<-read_tsv("./GPL6480-9577.txt", comment = "#")%>%
  select("ID","GENE_SYMBOL","GENE")%>%
  filter(!(str_detect(GENE_SYMBOL,"\\/")))%>%
  rename(Gene_Symbol=GENE_SYMBOL)
vali_exp<-expr_norm %>%as.data.frame()%>%
  rownames_to_column("ID")%>%
  inner_join(gpl6480_anno[,1:2],by="ID")%>%
  dplyr::rename("gene"="Gene_Symbol")%>%
  dplyr::select(-"ID")%>%
  aggregate(.~ gene, data = ., max)%>%
  column_to_rownames("gene")
gse13213_exp<-vali_exp
gse13213_pd<-target
save(gse13213_exp,gse13213_norm,gse13213_pd,file="./01GSE13213_117Samples.rda")

rm(list=ls())
library(survival)
library(survminer)
library(tidyverse)
load("./01GSE13213_117samples.rda")
load("/Volumes/T7/Methylation/01Pheno_Probe.rda")
gene<-unique(pheno$gene)
target<-gse13213_pd
vali_exp<-gse13213_exp
table(target$Stage)
target$age_group[target$Age<65]<-"<65"
target$age_group[target$Age>=65]<-">=65"
target$stage_group[target$Stage%in%c("IA","IB","IIA","IIB")]<-"StageI_II"
target$stage_group[target$Stage%in%c("IIIA","IIIB")]<-"StageIII_IV"
colnames(target)
Clinic<-target[,c(24,6,8,18:20,23,15,16,17)]
Clinic$p53[Clinic$p53==names(table(Clinic$p53)[2])]<-NA
Clinic$Smoking<-ifelse(Clinic$Smoking==0,"No","Yes")
restab<-descrTable(data=Clinic,show.all = T)
export2csv(restab,file="./Validation_Result/00Clinic_feature_gse13213.csv")
## Validation according to driver genes status-----
library(tidyverse)
library(survival)
library(survminer)
library(survcomp)
load("/Users/huzixin/Data_Project/Methylation/00mIMg.rda")
mIMg<-result
load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
##01GSE11969--------
rm(list=ls())
load("./01GSE11969_90samples.rda")
vali_exp<-gse11969_exp
target<-gse11969_norm$targets
table(target$EGFR)
#Mut  Wt 
#32  58
# Mut  Wt 
# 45  72 
target$age_group[target$Age<65]<-"<65"
target$age_group[target$Age>=65]<-">=65"
table(target$Stage)
target$stage_group[target$Stage%in%c("IA","IB","IIA","IIB")]<-"StageI_II"
target$stage_group[target$Stage%in%c("IIIA","IIIB")]<-"StageIII_IV"


##02GSE13213--------
load("./01GSE13213_117samples.rda")
vali_exp<-gse13213_exp
target<-gse13213_pd

##common codes-----
gene<-unique(mIMg$gene)
Clinic<-target%>%
  filter(!is.na(time)&!(time==0)) %>% 
  filter(Kras=="Wt")
Surv_ob<-Surv(time=Clinic$time,event=Clinic$status) 
gene<-str_replace_all(gene,"-","_")
rownames(vali_exp)<-str_replace_all(rownames(vali_exp),"-","_")
exp<-vali_exp[,Clinic$Sample]%>%
  dplyr::filter(rownames(vali_exp)%in%gene)%>%
  t()%>%as.data.frame()
result<-data.frame(matrix(NA,ncol=5,nrow=length(colnames(exp))))
gene<-colnames(exp)
for(i in 1:length(colnames(exp))){
  dat<-cbind(t(vali_exp[gene[i],Clinic$Sample]),
             Clinic[,c("age_group","Sex","stage_group","status","time")]) %>% 
    as.data.frame() 
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  dat[,1]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(colnames(dat)[1:4],collapse="+"))) 
  GCox <- coxph(FML,data=dat) 
  GSum <- summary(GCox) 
  HR <- GSum$coefficients[,2][1]
  PValue <- GSum$coefficients[,5][1]
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result[i,]<-c(gene[i],HR,names(HR),CI,PValue)
}
colnames(result)<-c("gene","HR","group","95%CI","P.val")
result<-result %>% 
  arrange(P.val)
result$P.val<-as.numeric(result$P.val)
result$HR<-as.numeric(result$HR)
result_gpl750<-result
result_gpl6884<-result
result_gpl11969<-result
result_gpl13213<-result
save(result_gpl11969,result_gpl13213,file="./Validation_Result/egfr_mut_gpl_multicox.rda")
save(result_gpl11969,result_gpl13213,file="./Validation_Result/tp53_wt_gpl_multicox.rda")
save(result_gpl11969,result_gpl13213,file="./Validation_Result/kras_wt_gpl_multicox.rda")

## EGFR/SRGN Figure 7 F-G------
dat<-vali_exp['SRGN',target$Sample] %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status","EGFR")])

dat1<-dat[dat$EGFR=="Wt",]
dat2<-dat[dat$EGFR=="Mut",]
res.cut <- surv_cutpoint(dat2, time = "time", event = "status",
                         variables = 'SRGN' )
quantile(dat2$SRGN)
summary(res.cut)
res.cat <- surv_categorize(res.cut)
dat2$exp_group<-res.cat$SRGN

Surv_ob<-Surv(dat2$time/30, dat2$status)
FML<-as.formula(paste0("Surv_ob","~",'exp_group'))
fit <- survfit(FML, data = dat2)
Fit<-coxph(FML, data = dat2) 
survdiff(FML, data = dat2)
modelSum <- summary(Fit) ;
modelSum$conf.int
# gse11969
#wt
# exp(coef) exp(-coef) lower .95 upper .95
# exp_grouplow 0.3769177   2.653099 0.1290386  1.100965
#0.38(0.13,1.10)
#             exp(coef) exp(-coef) lower .95 upper .95
#exp_grouplow  5.397676  0.1852649  1.882383  15.47767
#5.40 (1.88 - 15.48)
#gse13213
#exp(coef) exp(-coef) lower .95 upper .95
#exp_grouplow 0.3210783   3.114505 0.1365668 0.7548779
#0.32(0.14-0.75) 

cindex <- concordance.index(predict(Fit),surv.time = dat2$time, surv.event = dat2$status ,method = "noether") 
cindex$c.index
cindex$lower
cindex$upper
paste0(round(cindex$c.index,2)," (",round(cindex$lower,2)," - ",round(cindex$upper,2),")")
#gse11969
#wt
# 0.74(0.55,0.94)
#mut
# 0.86(0.74,0.98)
#gse13213
#wt
#0.75 (0.58 - 0.91)
#mut
#0.8 (0.65 - 0.95)
dir.create("./Validation_Result/EGFR_SRGN/")
pdf("./Validation_Result/EGFR_SRGN/gse13213_mut.pdf",height=6,width = 6)
ggsurvplot(survfit(Surv(dat2$time/30, dat2$status) ~ exp_group,data=dat2), data=dat2,
           conf.int=F, pval=T,pval.method=T,risk.table=T,
           ggtheme=theme_classic(),
           surv.median.line="hv",risk.table.height=0.25,
           risk.table.y.text.col=T,
           risk.table.y.text=F,
           ncensor.plot=F,font.x=c(12,'plain','black'),
           font.y=c(12,'plain','black'),
           font.tickslab=c(12,'plain','black'),
           xlab="Time in months",palette="Set1",
           legend.title = 'SRGN',
           legend.labs = c("High", "Low"),legend=c(0.8,0.8),
           xlim=c(0,180),break.time.by=50)
dev.off()
## EGFR/SRGN ROC Figure 7 J-K-----
dat1<-vali_exp['SRGN',target$Sample] %>% t() %>% as.data.frame()
dat<-dat1 %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status","EGFR","stage_group","age_group","Sex")]) %>% 
  filter(EGFR=="Wt")
variable<-'SRGN'

factor<-variable
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor)
res.cat <- surv_categorize(res.cut)
dat[,factor]<-res.cat[,3]

dat<-dat %>% filter(!(is.na(stage_group)) & !(is.na(age_group)) & !(is.na(Sex)))

Surv_ob<-Surv(dat$time/30, dat$status)

FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','stage_group',' + ','age_group',' + ','Sex')))
cox.timp<- coxph(FML,data = dat)
lpFit <- cox.timp$linear.predictors
roc.fit<-timeROC(T =dat$time/30,
                 delta=dat$status,
                 marker = lpFit, 
                 cause = 1, 
                 weighting = "marginal", 
                 times = c(12*c(1,3,5,7,9)), 
                 ROC = T,
                 iid = TRUE) 

pdf(paste0("./Validation_Result/gpl570/Multi_Diagnosis/SRGN_gpl11969",factor,"_ROC.pdf"),width=4,height=4)
plot(roc.fit,time=12*5,col = "#ED0000FF",add = FALSE,title = F)
plot(roc.fit,time=12*7,col = "#42B540FF",add = TRUE)
plot(roc.fit,time=12*9,col = "#925E9FFF",add = TRUE)

auc<-round(roc.fit[["AUC"]],2)
legend("bottomright",
       legend = c(
         paste0('AUC at 5 years',' ',auc[3]),
         paste0('AUC at 7 years',' ',auc[4]),
         paste0('AUC at 9 years',' ',auc[5])),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,seg.len=1,cex=0.9,
       col = c("#ED0000FF","#42B540FF","#925E9FFF"),
       bty = 'n')
title(main=factor,)
dev.off()

## EGFR/SRGN expression Figure 7 N-O-------
library(survival)
library(survminer)
rm(list=ls())
load("./01GSE13213_117samples.rda")
vali_exp<-gse13213_exp
target<-gse13213_pd
colnames(target)[15]<-"EGFR_status"
dat<-vali_exp[c('SRGN','EGFR'),target$Sample] %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","EGFR_status","time","status")],by="Sample")
dat1<-dat[dat$EGFR_status=="Wt",]
dat2<-dat[dat$EGFR_status=="Mut",]
res.cut1 <- surv_cutpoint(dat1, time = "time", event = "status",
                          variables = 'SRGN' )
res.cut2 <- surv_cutpoint(dat2, time = "time", event = "status",
                          variables = 'SRGN' )

summary(res.cut2)
res.cat1 <- surv_categorize(res.cut1)
res.cat2 <- surv_categorize(res.cut2)

dat1$exp_group<-res.cat1$SRGN
dat2$exp_group<-res.cat2$SRGN
dat<-rbind(dat1,dat2) %>% 
  as.data.frame() %>% 
  mutate(group=paste0(EGFR_status,"_",exp_group))
table(dat$group)
dat$group[dat$group=="Wt_low"]<-"group1"
dat$group[dat$group=="Wt_high"]<-"group2"
dat$group[dat$group=="Mut_low"] <-"group3"
dat$group[dat$group=="Mut_high"]<-"group4"

with(dat,tapply(SRGN,group,shapiro.test))
wilcox.test(SRGN~EGFR_status,data=dat[dat$exp_group=="low",],exact=F,alternative='two.sided')
#W = 393.5, p-value =0.0005863
wilcox.test(SRGN~EGFR_status,data=dat[dat$exp_group=="high",],exact=F,alternative='two.sided')
#W = 4, p-value = 4.64e-05
wilcox.test(SRGN~EGFR_status,data=dat,exact=F,alternative='two.sided')

library(RColorBrewer)
color<-brewer.pal(4,"Set2")
names(color)<-c("group1","group2","group3","group4")
colnames(dat_im)
my_comparisons<-list(c("group1,group3"),c("group2","group4"))
pdf("./Validation_Result/EGFR_SRGN/SRGN_revised_byegfr.pdf",height=3,width=4)
ggplot(dat ,aes(EGFR_status,SRGN,fill=group))+
  geom_violin()+
  theme_classic()+
  facet_wrap("exp_group",ncol=2)+
  #ylim(9.8,10.6)+
  labs(x = "", y = "Expression of SRGN")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=color)+
  #scale_y_continuous()+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12),
        strip.background=element_blank(),
        strip.text = element_text(colour = "black",size=12))
#stat_compare_means(label="p.signif",method="t.test")+  
#stat_compare_means()
dev.off()

cor.test(dat$EGFR[dat$exp_group=="high"],dat$SRGN[dat$exp_group=="high"],method="spearman")

cor.test(dat$EGFR[dat$exp_group=="low"],dat$SRGN[dat$exp_group=="low"],method="spearman")
dat$SRGN

high
cor:  rho -0.4638282 p-value = 1.44e-05
low 
cor:  rho 0.04860008  p-value = 0.3597

pdf("./Validation_Result/EGFR_SRGN/CoreelationofSRGNEGFR_revised.pdf",height=5,width=4)
ggplot(dat,aes(x=EGFR,y=SRGN))+
  geom_jitter(size = 2,shape = 21, color = 'black',aes(fill = group))+
  geom_smooth(method=lm,formula= y ~ x)+ 
  labs(x = "Expression of EGFR", y = "Expression of SRGN")+
  theme_bw()+
  theme(legend.position="none")+
  scale_color_manual(values=color)+
  facet_wrap("exp_group",ncol=1,scales="free_y")+
  theme(axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=12))+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12),
        strip.background=element_blank(),
        strip.text = element_text(colour = "black",size=12))
dev.off()
