## Explore LCK CD247 PSTPIP1 cg09032544 cg07786657 cg11683242 cg26227523-----
## correlation Figure 5 A----
library(ggExtra)
library(tidyverse)
load("./02.TPM&Targets.Rda")
target<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
exp0<-log2(TPM_sym+1)
rm(TPM,TPM_sym,Targets)
load("./04CHAMPline_import.Rda")
beta0<-beta_norm
rm(beta_norm,Myload)
dat1<-exp0[c("LCK","CD247","PSTPIP1"),target$Sample] %>% t() %>% as.data.frame()
dat2<-beta0[c('cg09032544', 'cg07786657', 'cg11683242', 'cg26227523'),target$Sample] %>% t() %>% as.data.frame()
dat<-cbind(dat1,dat2)
dat<-dat[,c("cg09032544","cg07786657","CD247","cg11683242","LCK","cg26227523","PSTPIP1")]
load("02Cluster_IM.Rda")
CYT<-exp0[c("GZMA","PRF1"),]%>%
  t()%>%as.data.frame()%>%
  mutate(CYT=1/2*(log(PRF1)+log(GZMA)))
rownames(res_gsva)
Imm<-t(res_gsva)[target$Sample,c(2,3,4,19,28,30,21,22,15,16)]
Imm<-cbind(Imm[target$Sample,],
           CYT=CYT[target$Sample,"CYT"],
           PDCD1=t(exp0["PDCD1",target$Sample]))
rm(res_gsva,IPS,Type,CYT,IM,Targets)

library(corrplot)
corr<-WGCNA::cor(dat,Imm, use = "p",method="spearman")
nSamples<-nrow(dat)
p.mat<-WGCNA::corPvalueStudent(corr, nSamples)
dir.create("./fig5_mIMg/")
pdf("./fig5_mIMg/methy_Im_Corr.pdf",width=10,height=10)
corrplot(corr,p.mat=p.mat,sig.level=0.05,addshade = "all",
         addCoef.col = 'grey',col=pal_gsea("default", n = 100, reverse = F)(100),
         method = "circle",number.cex=1.2,tl.col="black",tl.cex=1)
dev.off()
dat<-as.data.frame(dat) %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","Cluster")],by="Sample")
cor.test(dat$cg26227523,dat$PSTPIP1,method="spearman")

## Figure5 B------
par(mfrow=c(1,1))
pdf("./fig5_mIMg/PSTPIP1_cg26227523.pdf",width=4,height=4)
p<-ggplot(dat,aes(x=cg26227523,y=PSTPIP1))+
  geom_point(size = 2,shape = 1, fill = 'gray',aes(color=Cluster))+
  geom_smooth(method=lm,formula= y ~ x)+
  theme_bw()+
  #theme(legend.position="none")+
  scale_color_manual(values=c('#48466d','#3d84a8','#46cdcf'))+
  annotate("text", label = paste0("Spearman correlation","\n","r = -0.64","\n","p < 2.2e-16"),
           x = 0.8, y = 5.5, size = 4, colour = "#696969",family = "serif", fontface = "italic")+
  theme(axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=12))
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)
dev.off()
legend<-get_legend(p)
pdf("./fig5_mIMg/legend.pdf",width=2,height=2)
plot_grid(legend)
dev.off()
## Figure6 A-G------
target<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
target0<-target
colnames(target)[c(9,11)]<-c("time","status")
target<-target %>% 
  filter(!is.na(time)&!(time==0))
dat<-cbind(dat1,dat2) %>% 
  select(c("cg09032544","cg07786657","CD247","cg11683242","LCK","cg26227523","PSTPIP1")) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status")],by="Sample")
colnames(dat)
variable<-colnames(dat)[2:8]
library(survminer)
library(survival)
library(survcomp)
result<-data.frame(matrix(NA,ncol=7))
colnames(result)<-c("p.val","HR","HR_lower","HR_upper","c-index","c-index_lower","c-index_upper")
for(i in 1:7){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor )
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  Surv_ob<-Surv(res.cat$time/30, res.cat$status)
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
write.table(result,"./fig5_mIMg/mIMg_survival.txt",row.names=T,col.names=T,quote=F,sep="\t")
factor<-variable[3]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time/30, res.cat$status)
FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)


pdf(paste0("./fig5_mIMg/",variable[7],"survival.pdf"),height=6,width = 6)
ggsurvplot(survfit(Surv(res.cat$time/30, res.cat$status)~PSTPIP1,data=res.cat), data=res.cat,
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


## Figure6 ROC H-N-----

library(timeROC)
timeROC_helper<-function(x){
  
  cox.timp<- coxph(Surv(res.cat$time/30,res.cat$status==1) ~ x,data = res.cat)
  lpFit <- cox.timp$linear.predictors
  roc.fit <-timeROC(T =res.cat$time,
                    delta=res.cat$status,
                    marker = lpFit, 
                    cause = 1, 
                    weighting = "marginal", 
                    times = c(12*c(1,3,5,7,9,15)), 
                    ROC = T,
                    iid = TRUE) 
}


par(mfrow=c(1,1))
par(mgp=c(1.6,1,0),mar=c(3,3,3,1)+0.1)
par(cex.axis=1,cex.lab=1,font.lab=1,font.axis=1,font.main=1)
factor<-variable[1]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )

summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time/30, res.cat$status)

FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)
cox.timp<- coxph(FML,data = res.cat)

roc.fit<-timeROC_helper(res.cat$cg09032544)
pdf(paste0("./fig5_mIMg/",factor,"_ROC.pdf"),width=4,height=4)
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

dev.off()
roc.fit$AUC

## Figure S6 A,F-I_ACAP1 ARHGAP30-------- 

dat1<-beta0[c("cg20009327","cg02358862","cg15698795"),target$Sample] %>% t() %>% as.data.frame()
dat2<-exp0[c("ARHGAP9","CORO1A","TBC1D10C"),target$Sample] %>% t() %>% as.data.frame()
cor.test(dat1$cg20009327,dat2$ARHGAP9,method="spearman")
cor.test(dat1$cg02358862,dat2$CORO1A,method="spearman")
cor.test(dat1$cg15698795,dat2$TBC1D10C,method="spearman")

dat1<-exp0[c("ACAP1","ARHGAP30"),target$Sample] %>% t() %>% as.data.frame()
dat2<-beta0[c("cg25671438","cg01774645"),target$Sample] %>% t() %>% as.data.frame()
dat<-cbind(dat1,dat2)
dat<-dat[,c("cg25671438","ACAP1","cg01774645","ARHGAP30")]
corr<-WGCNA::cor(dat,Imm, use = "p",method="spearman")
nSamples<-nrow(dat)
p.mat<-WGCNA::corPvalueStudent(corr, nSamples)
dir.create("./fig5_mIMg/")
pdf("./fig5_mIMg/Gap_methy_Im_Corr.pdf",width=10,height=10)
corrplot(corr,p.mat=p.mat,sig.level=0.05,addshade = "all",
         addCoef.col = 'grey',col=pal_gsea("default", n = 100, reverse = F)(100),
         method = "circle",number.cex=1.2,tl.col="black",tl.cex=1)
dev.off()


target0<-target
colnames(target)[c(9,11)]<-c("time","status")
target<-target %>% 
  filter(!is.na(time)&!(time==0))
dat<-dat %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status")],by="Sample")
colnames(dat)
variable<-colnames(dat)[2:6]
library(survminer)
library(survival)
library(survcomp)
result<-data.frame(matrix(NA,ncol=7))
colnames(result)<-c("p.val","HR","HR_lower","HR_upper","c-index","c-index_lower","c-index_upper")
for(i in 1:5){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor )
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  Surv_ob<-Surv(res.cat$time/30, res.cat$status)
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
write.table(result,"./fig5_mIMg/gap_mIMg_survival.txt",row.names=T,col.names=T,quote=F,sep="\t")
factor<-variable[5]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time/30, res.cat$status)
FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)
pdf(paste0("./fig5_mIMg/",variable[5],"survival.pdf"),height=6,width = 6)
ggsurvplot(survfit(Surv(res.cat$time/30, res.cat$status)~ARHGAP30,data=res.cat), data=res.cat,
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
par(mfrow=c(1,1))
par(mgp=c(1.6,1,0),mar=c(3,3,3,1)+0.1)
par(cex.axis=0.75,cex.lab=0.75,font.lab=1,font.axis=1,font.main=1)
factor<-variable[5]
res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                         variables = factor )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
Surv_ob<-Surv(res.cat$time/30, res.cat$status)
FML<-as.formula(paste0("Surv_ob","~",factor))
fit <- survfit(FML, data = res.cat)
roc.fit<-timeROC_helper(res.cat$ARHGAP30)
pdf(paste0("./fig5_mIMg/",factor,"_ROC.pdf"),width=4,height=4)
plot(roc.fit,time=12*5,col = "#00468BFF",add =FALSE,title = F)
plot(roc.fit,time=12*7,col = "#ED0000FF",add = TRUE)
plot(roc.fit,time=12*9,col = "#42B540FF",add = TRUE)
auc<-round(roc.fit[["AUC"]],2)
legend("bottomright",
       legend = c(paste0('AUC at 5 year',' ',auc[3]), 
                  paste0('AUC at 7 years',' ',auc[4]),
                  paste0('AUC at 9 years',' ',auc[5])),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,seg.len=1,cex=0.8,
       col = c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF"),
       bty = 'n')
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



## Figure S5 A-C Prognostic model Diagnosis DCA  -----
library(ggExtra)
library(survminer)
library(tidyverse)
library(survMisc)
library(foreign)
library(Hmisc)
library(rms)
library(survival)
library(timeROC)
source("/Users/huzixin/Data_Project/Rscript/DCA/stdca.R")
load("./02.TPM&Targets.Rda")
target<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
exp0<-log2(TPM_sym+1)
rm(TPM,TPM_sym,Targets)
load("./04CHAMPline_import.Rda")
beta0<-beta_norm
rm(beta_norm,Myload)
dat1<-exp0[c("LCK","CD247","PSTPIP1"),target$Sample] %>% t() %>% as.data.frame()
dat2<-beta0[c('cg09032544', 'cg07786657', 'cg11683242', 'cg26227523'),target$Sample] %>% t() %>% as.data.frame()
dat<-cbind(dat1,dat2)
target<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
target0<-target
colnames(target)[c(9,11)]<-c("time","status")
target<-target %>% 
  filter(!is.na(time)&!(time==0))
dat<-cbind(dat1,dat2) %>% 
  select(c("cg09032544","cg07786657","CD247","cg11683242","LCK","cg26227523","PSTPIP1")) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status","stage_group","age_group","Gender")],by="Sample")
colnames(dat)
variable<-colnames(dat)[2:8]
for(i in 1:length(variable)){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor )
  res.cat <- surv_categorize(res.cut)
  dat[,factor]<-res.cat[,3]
}
dat<-dat %>% filter(!(is.na(stage_group)) & !(is.na(age_group)) & !(is.na(Gender)))

Surv_ob<-Surv(dat$time/30, dat$status)
fml_c<-as.formula(paste0("Surv_ob","~",paste0('stage_group',' + ','age_group',' + ','Gender')))
coxmod_c<-coxph(fml_c,data=dat)

dir.create("./Multi_Diagnosis/DCA/")
par(cex.axis=1,cex.lab=1,font.lab=1,font.axis=1,font.main=1)
for(i in 1:length(variable)){
  factor<-variable[i]
  FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','stage_group',' + ','age_group',' + ','Gender')))
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
  
  pdf(paste0("./Multi_Diagnosis/ROC_TCGA/",factor,"_ROC.pdf"),width=4,height=4)
  #plot(roc.fit,time=12,col = "#00468BFF",add =FALSE,title = F)
  #plot(roc.fit,time=12*3,col = "#ED0000FF",add = TRUE)
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


variate0<-c("CD247","LCK","PSTPIP1")
for(i in 1:3){
  factor<-variate0[i]
  fml_c<-as.formula(paste0("Surv_ob","~",factor))
  coxmod_c<-coxph(fml_c,data=dat)
  dat$signature_c<-c(1-(summary(survfit(coxmod_c,newdata=dat),times=12*3)$surv))
  km_c<-stdca(data=dat,outcome="status",ttoutcome="time",timepoint=365*3,
              predictors="signature_c",loess.span = 0.5,smooth=T)
  FML<-as.formula(paste0("Surv_ob","~",paste0(factor,' + ','stage_group',' + ','age_group',' + ','Gender')))
  coxmod<-coxph(FML,data=dat)
  dat$signature<-c(1-(summary(survfit(coxmod,newdata=dat),times=12*3)$surv))
  km<-stdca(data=dat,
            outcome="status",ttoutcome="time",
            timepoint= 365*3,
            predictors="signature",
            loess.span = 0.5,smooth=T)
  x_max<-round(max(km$net.benefit$threshold),1)
  x_min<-round(min(km$net.benefit$threshold),1)
  y_max<-round(max(km$net.benefit$all),1)
  y_min<--0.05
  pdf(paste0("./Multi_Diagnosis/DCA/",factor,"_3years","_DCA.pdf"),width=4,height=4)
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