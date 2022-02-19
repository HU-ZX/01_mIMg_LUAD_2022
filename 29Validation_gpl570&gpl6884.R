library(GEOquery)
library(tidyverse)
library(sva)
library(GSVA)
library(ggpubr)
library(ggsci)
library(survival)
library(survminer)
library(rms)
library(survMisc)
library(foreign)
library(Hmisc)
library(rms)
library(survival)
library(survcomp)
library(timeROC)

###gpl570------
load("./01GSE14814_68samples.rda")
load("./01GSE50081_112samples.rda")
load("./01GSE37745_92samples.rda")
exp_37745<-exprs(gse37745_rma0)
exp_50081<-exprs(gse50081_rma0)
exp_14814<-exprs(gse14814_rma)%>%as.data.frame()%>%rownames_to_column("ID")
identical(rownames(exp_37745),rownames(exp_50081))
cb<-cbind(exp_37745,exp_50081)%>%as.data.frame()%>%
  rownames_to_column("ID")%>%
  inner_join(exp_14814,by="ID")%>%
  column_to_rownames("ID")%>%
  as.matrix()
cb2<-cbind(exp_37745,exp_50081)

target_37745<-pData(gse37745_rma0)%>%
  dplyr::select(c("Sample","status","time"))
target_50081<-pData(gse50081_rma0)%>%
  dplyr::select(c("Sample","status","time"))
target_14814<-pData(gse14814_rma)%>%
  dplyr::select(c("Sample","status","time"))
target_37745$GSE<-rep("GSE37745",nrow(target_37745))
target_50081$GSE<-rep("GSE50081",nrow(target_50081))
target_14814$GSE<-rep("GSE14814",nrow(target_14814))
target<-bind_rows(target_37745,target_50081,target_14814)
batch<-target$GSE
target2<-bind_rows(target_37745,target_50081)
batch2<-target2$GSE
library(sva)
combat_exp <- ComBat(dat = cb,
                     batch = batch,
                     # mod = mod,
                     par.prior = TRUE,
                     prior.plots = FALSE)
combat_exp2 <- ComBat(dat = cb2,
                      batch = batch2,
                      # mod = mod,
                      par.prior = TRUE,
                      prior.plots = FALSE)
gpl570_anno<-read_tsv("./GPL570-55999.txt", comment = "#")%>%
  select("ID","Gene Symbol","ENTREZ_GENE_ID")%>%
  filter(!(str_detect(`Gene Symbol`,"\\/")))%>%
  rename("gene"="Gene Symbol")

vali_exp<-combat_exp%>%as.data.frame()%>%
  rownames_to_column("ID")%>%
  inner_join(gpl570_anno[,1:2],by="ID")%>%
  dplyr::select(-"ID")%>%
  aggregate(.~ gene, data = ., max)%>%
  column_to_rownames("gene")
vali_exp2<-combat_exp2%>%as.data.frame()%>%
  rownames_to_column("ID")%>%
  inner_join(gpl570_anno[,1:2],by="ID")%>%
  dplyr::select(-"ID")%>%
  aggregate(.~ gene, data = ., max)%>%
  column_to_rownames("gene")
save(vali_exp,vali_exp2,target,target2,file="./gpl570_cbdata.rda")
library(compareGroups)
colnames(target)
Clinic<-target[,c("age","gender","Stage","T","N","M","tobacco_history","GSE")]
restab<-descrTable(data=Clinic,show.all = T)
export2csv(restab,file="./Validation_Result/gpl570/00Clinic_feature_gpl570.csv")
###gpl6884-----
rm(list=ls())
library(tidyverse)
library(survival)
library(survminer)
load("./01GSE41271_183samples.rda")
load("./01GSE42127_133samples.rda")
load("/Volumes/T7/hzx/data_immune_non_diploidy/Validation/GPL/GPL_annobam.Rda")
GPL6884_table<-GPL6884_entrz
rm(GPL6480_bam,GPL6480_entrz,
   GPL570_bam,GPL570_entrz,
   GPL6884_bam,GPL6884_entrz)
target_41271 <-  gse41271_pd  
target_41271 <- gse41271_pd[colnames(gse41271_expr),]
target_41271$GSE<-rep("GSE41271",nrow(target_41271))
colnames(gse41271_expr)<-target_41271 $Sample

target_42127 <-  gse42127_pd  
target_42127$GSE<-rep("GSE42127",nrow(target_42127))
target<-bind_rows(target_41271,target_42127)
exp_41271<-gse41271_expr%>%as.data.frame()%>%rownames_to_column("probe_id")
exp_42127<-gse42127_exp%>%as.data.frame()%>%rownames_to_column("probe_id")
vali_exp<-exp_41271%>%inner_join(exp_42127,by="probe_id")%>%
  inner_join(GPL6884_table[,c(1,3)],by="probe_id")%>%
  dplyr::rename("gene"="symbol")%>%
  dplyr::select(-"probe_id")%>%
  aggregate(.~ gene, data = ., max)%>%
  column_to_rownames("gene")
save(vali_exp,target,file="./gpl6884_cbdata.rda")
colnames(target)
Clinic<-target[,c("age","gender","Stage","tobacco_history","GSE")]
restab<-descrTable(data=Clinic,show.all = T)
export2csv(restab,file="./Validation_Result/gpl6884/00Clinic_feature_gpl6884.csv")


## common code FigureS4-5-----
## Validation of prognostic mIMg -----
rm(list=ls())
load("./gpl6884_cbdata.rda")
load("/Users/huzixin/Data_Project/Methylation/00mIMg.rda")
gene<-unique(result$gene)
gene<-str_replace_all(gene,"-","_")
rownames(vali_exp)<-str_replace_all(rownames(vali_exp),"-","_")
Clinic<-target%>%
  filter(!is.na(time)&!(time==0))
Surv_ob<-Surv(time=Clinic$time,event=Clinic$status) 
exp<-vali_exp[,Clinic$Sample]%>%
  dplyr::filter(rownames(vali_exp)%in%gene)%>%
  t()%>%as.data.frame()
colnames(exp)
gene<-setdiff(gene,colnames(exp))
result<-data.frame(matrix(NA,ncol=5,nrow=length(colnames(exp))))
gene<-colnames(exp)
for(i in 1:length(colnames(exp))){
  dat<-cbind(t(vali_exp[gene[i],Clinic$Sample]),
             Clinic[,c("age","gender","Stage","status","time")]) %>% 
    as.data.frame() 
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  dat[,1]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(colnames(dat)[1:4],collapse="+"))) 
  GCox <- coxph(FML,data=dat) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2][1],2)
  PValue <- round(GSum$coefficients[,5][1],3)
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result[i,]<-c(gene[i],HR,names(HR),CI,PValue)
}
colnames(result)<-c("gene","HR","group","95%CI","P.val")
result$HR<-as.numeric(result$HR)
result$P.val<-as.numeric(result$P.val)
result0<-result
result_gpl570<-rbind(result0,result)
result_gpl6884<-result
save(result_gpl570,result_gpl6884,file="./Validation_Result/gpl_multicox.rda")

Unicox_result<-data.frame(matrix(NA,ncol=4,nrow=length(colnames(exp))))
gene<-colnames(exp)
for(i in 1:length(colnames(exp))){
  dat<-cbind(t(vali_exp[gene[i],Clinic$Sample]),
             Clinic[,c("status","time")]) %>% 
    as.data.frame() 
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  dat[,1]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~',gene[i],collapse="+")) 
  GCox <- coxph(FML,data=dat) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") 
  Unicox_result[i,]<-c(gene[i],HR,CI,PValue)
}
colnames(Unicox_result)<-c("gene","HR","95%CI","P.val")
Unicox_result$HR<-as.numeric(Unicox_result$HR) 
Unicox_result$P.val<-as.numeric(Unicox_result$P.val) 

Unicox_result0<-Unicox_result
Unicox_result_gpl570<-rbind(Unicox_result0,Unicox_result[2,])
Unicox_result_gpl6884<-Unicox_result
save(Unicox_result_gpl570,Unicox_result_gpl6884,file="./Validation_Result/gpl_unicox.rda")

## survival curves------
dat1<-vali_exp[c("LCK","CD247","PSTPIP1"),target$Sample] %>% t() %>% as.data.frame()
dat1<-vali_exp[c("ACAP1","ARHGAP30"),target$Sample] %>% t() %>% as.data.frame()

dat<-dat1 %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c('time','status','Sample')],by="Sample")%>% 
  filter(!is.na(time)&!(time==0))
variable<-colnames(dat)[2:3]
library(survminer)
library(survival)
library(survcomp)
result<-data.frame(matrix(NA,ncol=7))
colnames(result)<-c("p.val","HR","HR_lower","HR_upper","c-index","c-index_lower","c-index_upper")
for(i in 1:2){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor )
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  Surv_ob<-Surv(res.cat$time, res.cat$status)
  FML<-as.formula(paste0("Surv_ob","~",factor))
  fit <- survfit(FML, data = res.cat)
  Fit<-coxph(FML, data = res.cat) 
  modelSum <- summary(Fit) ;
  HR<-round(modelSum$conf.int[1],2)
  HR_lower<-round(modelSum$conf.int[3],2)
  HR_upper<-round(modelSum$conf.int[4],2)
  cindex <- concordance.index(predict(Fit),surv.time = res.cat$time, surv.event = res.cat$status ,method = "noether") 
  result[i,]<-c(round(modelSum$sctest[3],3),HR,HR_lower,HR_upper,round(cindex$c.index,2),round(cindex$lower,2), round(cindex$upper,2))
  
}

rownames(result)<-variable
write.table(result,"./Validation_Result/gpl6884/survival/vali_mIMg_survival.txt",row.names=T,col.names=T,quote=F,sep="\t")

factor<-variable[2]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time, res.cat$status)
FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)
pdf(paste0("./Validation_Result/gpl6884/survival/",variable[2],"survival.pdf"),height=6,width = 6)
ggsurvplot(survfit(Surv(res.cat$time, res.cat$status)~ARHGAP30,data=res.cat), data=res.cat,
           conf.int=F, pval=T,pval.method=T,risk.table=T, 
           ggtheme=theme_classic(),
           surv.median.line="hv",risk.table.height=0.25,
           risk.table.y.text.col=T,
           risk.table.y.text=F,
           ncensor.plot=F,font.x=c(12,'plain','black'),
           font.y=c(12,'plain','black'),
           font.tickslab=c(12,'plain','black'),
           xlab="Time in months",palette="Set1",
           legend.title = factor,
           legend.labs = c("Hyper", "Hypo"),legend=c(0.8,0.8),
           xlim=c(0,180),break.time.by=50)
dev.off()
## ROC curves------
library(timeROC)
timeROC_helper<-function(x){
  
  cox.timp<- coxph(Surv(res.cat$time/30,res.cat$status==1) ~ x,data = res.cat)
  lpFit <- cox.timp$linear.predictors
  roc.fit <-timeROC(T =res.cat$time,
                    delta=res.cat$status,
                    marker = lpFit, 
                    cause = 1, 
                    weighting = "marginal", 
                    times = c(12*c(1,3,5,7)), 
                    ROC = T,
                    iid = TRUE) 
}


par(mfrow=c(1,1))
par(mgp=c(1.6,1,0),mar=c(3,3,3,1)+0.1)
par(cex.axis=0.75,cex.lab=0.75,font.lab=1,font.axis=1,font.main=1)
factor<-variable[1]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time, res.cat$status)
FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)
roc.fit<-timeROC_helper(res.cat$ARHGAP30)
pdf(paste0("./Validation_Result/gpl6884/survival/",factor,"_ROC.pdf"),width=4,height=4)
plot(roc.fit,time=12,col = "#00468BFF",add =FALSE,title = F)
plot(roc.fit,time=12*3,col = "#ED0000FF",add = TRUE)
plot(roc.fit,time=12*5,col = "#42B540FF",add = TRUE)
plot(roc.fit,time=12*7,col = "#0099B4FF",add = TRUE)

auc<-round(roc.fit[["AUC"]],2)
legend("bottomright",
       legend = c(paste0('AUC at 1 year',' ',auc[1]), 
                  paste0('AUC at 3 years',' ',auc[2]),
                  paste0('AUC at 5 years',' ',auc[3]),
                  paste0('AUC at 7 years',' ',auc[4])),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,seg.len=1,cex=0.8,
       col = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"),
       bty = 'n')
title(factor)
dev.off()
roc.fit$AUC



## Model Diagnosis DCA curves------
rm(list=ls())
source("/Users/huzixin/Data_Project/Rscript/DCA/stdca.R")
load("./gpl6884_cbdata.rda")
dat1<-vali_exp[c("LCK","CD247","PSTPIP1"),target$Sample] %>% t() %>% as.data.frame()
dat<-dat1 %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c('time','status','Sample','age','Stage','gender')],by="Sample") 
variable<-colnames(dat)[2:4]

for(i in 1:length(variable)){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor )
  res.cat <- surv_categorize(res.cut)
  dat[,factor]<-res.cat[,3]
}

dat<-dat %>% filter(!(is.na(Stage)) & !(is.na(age)) & !(is.na(gender)))
Surv_ob<-Surv(dat$time, dat$status)
dir.create("./Validation_Result/gpl6884/Multi_Diagnosis/DCA")
par(cex.axis=1,cex.lab=1,font.lab=1,font.axis=1,font.main=1)
for(i in 1:length(variable)){
  factor<-variable[i]
  FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','Stage',' + ','age',' + ','gender')))
  cox.timp<- coxph(FML,data = dat)
  lpFit <- cox.timp$linear.predictors
  roc.fit<-timeROC(T =dat$time,
                   delta=dat$status,
                   marker = lpFit, 
                   cause = 1, 
                   weighting = "marginal", 
                   times = c(12*c(1,3,5,7,9)), 
                   ROC = T,
                   iid = TRUE) 
  
  pdf(paste0("./Validation_Result/gpl6884/Multi_Diagnosis/ROC/",factor,"_ROC.pdf"),width=4,height=4)
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
}


for(i in 1:length(variable)){
  factor<-variable[i]
  FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','Stage',' + ','age',' + ','gender')))
  coxmod<-coxph(FML,data=dat)
  dat$signature<-c(1-(summary(survfit(coxmod,newdata=dat),times=12*3)$surv))
  km<-stdca(data=dat,
            outcome="status",ttoutcome="time",
            timepoint= 12*3,
            predictors="signature",
            loess.span = 0.5,smooth=T)
  x_max<-round(max(km$net.benefit$threshold),1)
  x_min<-round(min(km$net.benefit$threshold),1)
  y_max<-round(max(km$net.benefit$all),1)
  y_min<-0
  #dat$signature_c<-c(1-(summary(survfit(coxmod_c,newdata=dat),times=12*(2*j-1))$surv))
  # km_c<-stdca(data=dat,outcome="status",ttoutcome="time",timepoint=365*(2*j-1),
  #             predictors="signature_c",loess.span = 0.5,smooth=T)
  pdf(paste0("./Validation_Result/gpl6884/Multi_Diagnosis/DCA/",factor,"_3years","_DCA.pdf"),width=4,height=4)
  plot(km$net.benefit.threshold, km$net.benefit.none, 
       type = "l", lwd=2, xlim=c(x_min,x_max), ylim=c(y_min, y_max), 
       xlab = "Threshold Probability", ylab = "Net Benefit")
  lines(km$net.benefit$threshold, km$net.benefit$all, type="l", col=3, lwd=2)
  lines(km$net.benefit$threshold, km$net.benefit$signature, type="l", col=10 , lwd = 2)
  #lines(km_c$net.benefit$threshold, km_c$net.benefit$signature_c, type="l", col = 12, lty=2, lwd = 2)
  abline(h = 0, lwd =2,col  ="black")
  legend("topright", cex=1, legend=c("None", "All", factor), col=c(17, 3, 10), lwd=c(2, 2, 2), lty=c(1, 1, 1), bty = "n")
  dev.off()
}

for(i in 1:3){
  factor<-variable[i]
  fml_c<-as.formula(paste0("Surv_ob","~",factor))
  coxmod_c<-coxph(fml_c,data=dat)
  dat$signature_c<-c(1-(summary(survfit(coxmod_c,newdata=dat),times=12*3)$surv))
  km_c<-stdca(data=dat,outcome="status",ttoutcome="time",timepoint=12*3,
              predictors="signature_c",loess.span = 0.5,smooth=T)
  FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','Stage',' + ','age',' + ','gender')))
  coxmod<-coxph(FML,data=dat)
  dat$signature<-c(1-(summary(survfit(coxmod,newdata=dat),times=12*3)$surv))
  km<-stdca(data=dat,
            outcome="status",ttoutcome="time",
            timepoint= 12*3,
            predictors="signature",
            loess.span = 0.5,smooth=T)
  x_max<-round(max(km$net.benefit$threshold),1)
  x_min<-round(min(km$net.benefit$threshold),1)
  y_max<-round(max(km$net.benefit$all),1)
  y_min<--0.05
  pdf(paste0("./Validation_Result/gpl6884/Multi_Diagnosis/DCA/",factor,"_3years","_DCA.pdf"),width=4,height=4)
  plot(km$net.benefit.threshold, km$net.benefit.none, 
       type = "l", lwd=2, xlim=c(x_min,x_max), ylim=c(y_min, y_max), 
       xlab = "Threshold Probability", ylab = "Net Benefit")
  abline(h = 0, lwd =2,col  ="black")
  lines(km$net.benefit$threshold, km$net.benefit$all, type="l", col=3, lwd=2,lty=2)
  lines(km$net.benefit$threshold, km$net.benefit$signature, type="l", col=10 , lwd = 2)
  lines(km_c$net.benefit$threshold, km_c$net.benefit$signature_c, type="l", col = 12,  lwd = 2)
  legend("topright", cex=1, legend=c("None", "All", "Multivariate",factor), col=c(17, 3, 10,12), lwd=c(2, 2, 2,2), lty=c(1,2, 1, 1), bty = "n")
  dev.off()
}



