##EGFR------
rm(list=ls())
library(ggExtra)
library(tidyverse)
library(survminer)
load("./02.TPM&Targets.Rda")
target0<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
colnames(target0)[c(9,11)]<-c("time","status")
target<-target0%>% 
  filter(!(is.na(egfr_status))) %>% 
  filter(!is.na(time)&!(time==0))
exp0<-log2(TPM_sym+1)
rm(TPM,TPM_sym,Targets)
load("./04CHAMPline_import.Rda")
beta0<-beta_norm
rm(beta_norm,Myload)
dat1<-exp0['SRGN',target$Sample] %>% t() %>% as.data.frame()
dat2<-data.frame(cg02851793= beta0['cg02851793',target$Sample])
dat<-cbind(dat1,dat2) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c('Sample','time','status','egfr_status','Cluster')],by="Sample")



## Figure7 D-E Logrank -----
dat<-dat %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status","Cluster","stage_group","age_group","Gender")])

dat1<-dat[dat$egfr_status=="egfr_wt",]
dat2<-dat[dat$egfr_status=="egfr_mut",]
res.cut <- surv_cutpoint(dat2, time = "time", event = "status",
                         variables = 'cg02851793' )
res.cut <- surv_cutpoint(dat2, time = "time", event = "status",
                         variables = 'SRGN' )

summary(res.cut)
res.cat <- surv_categorize(res.cut)
dat2$methy_group<-res.cat$cg02851793
dat2$exp_group<-res.cat$SRGN

Surv_ob<-Surv(dat2$time/30, dat2$status)
FML<-as.formula(paste0("Surv_ob","~",'methy_group'))
FML<-as.formula(paste0("Surv_ob","~",'exp_group'))
fit <- survfit(FML, data = dat2)
Fit<-coxph(FML, data = dat2) 
survdiff(FML, data = dat2)
modelSum <- summary(Fit) ;
modelSum$conf.int
#                exp(coef) exp(-coef) lower .95 upper .95
# cg02851793low  1.964032  0.5091567  1.292271  2.984993
#          exp(coef) exp(-coef) lower .95 upper .95
# SRGNlow 0.5221187   1.915273 0.3484316 0.7823857

# mut
#          exp(coef) exp(-coef) lower .95 upper .95
#methy_grouplow 0.5405534   1.849956 0.2302534  1.269028
#exp_grouplow  3.397567  0.2943283  1.368043  8.437937

cindex <- concordance.index(predict(Fit),surv.time = dat2$time, surv.event = dat2$status ,method = "noether") 
cindex$c.index
cindex$lower
cindex$upper
#      0.65(0.55-0.76)
# SRGN 0.64(0.53-0.74)

# 0.68(0.47-0.89)
# 0.76(0.58-0.94)
pdf("./EGFR_SRGN/mut_methylation.pdf",height=6,width = 6)
ggsurvplot(survfit(Surv(dat2$time/30, dat2$status) ~ methy_group,data=dat2), data=dat2,
           conf.int=F, pval=T,pval.method=T,risk.table=T,
           ggtheme=theme_classic(),
           surv.median.line="hv",risk.table.height=0.25,
           risk.table.y.text.col=T,
           risk.table.y.text=F,
           ncensor.plot=F,font.x=c(12,'plain','black'),
           font.y=c(12,'plain','black'),
           font.tickslab=c(12,'plain','black'),
           xlab="Time in months",palette="Set1",
           legend.title = 'cg02851793',
           legend.labs = c("Hyper", "Hypo"),legend=c(0.8,0.8),
           xlim=c(0,180),break.time.by=50)
dev.off()

pdf("./EGFR_SRGN/mut_expression.pdf",height=6,width = 6)
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
           legend.labs = c("Hyper", "Hypo"),legend=c(0.8,0.8),
           xlim=c(0,180),break.time.by=50)
dev.off()




## Figure7 H-I EGFR/SRGN ROC -----
dat1<-exp0['SRGN',target$Sample] %>% t() %>% as.data.frame()
dat2<-data.frame(cg02851793= beta0['cg02851793',target$Sample])
dat<-cbind(dat1,dat2) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c("Sample","time","status","egfr_status","stage_group","age_group","Gender")]) %>% 
  filter(egfr_status=="egfr_wt")
variable<-colnames(dat)[2:3]
for(i in 1:2){
  factor<-variable[i]
  res.cut <- surv_cutpoint(dat, time = "time", event = "status",
                           variables = factor)
  res.cat <- surv_categorize(res.cut)
  dat[,factor]<-res.cat[,3]
}
dat<-dat %>% filter(!(is.na(stage_group)) & !(is.na(age_group)) & !(is.na(Gender)))

Surv_ob<-Surv(dat$time/30, dat$status)
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
  title(main=factor)
  dev.off()
  
}








## Figure8 D EGFR_Enrichment_Heatmap ------
library(msigdbr)
library(GSVA)
library(survminer)
library(ggExtra)
library(tidyverse)
library(limma)

load("./02.TPM&Targets.Rda")
target0<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster)
colnames(target0)[c(9,11)]<-c("time","status")
target<-target0%>% 
  filter(!(is.na(egfr_status))) %>% 
  filter(!is.na(time)&!(time==0))
exp0<-log2(TPM_sym+1)
rm(TPM,TPM_sym,Targets)
msigdb_GMTs <- "./"
msigdb <- "c5.go.v7.4.symbols.gmt.txt"
geneset <- getGmt(file.path(msigdb_GMTs, msigdb)) 
exp<-exp0[,target$Sample] %>% as.matrix()
es_all<- gsva(exp,geneset, method = "gsva", kcdf="Gaussian")
save(es_all,file="./SRGN_gsva.rda")
rm(list=ls())
load("./SRGN_gsva.rda")

dat1<-exp0['SRGN',target$Sample] %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c('Sample','time','status','egfr_status')],by="Sample") %>% 
  filter(egfr_status=="egfr_wt")
dat2<-exp0['SRGN',target$Sample] %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(target[,c('Sample','time','status','egfr_status')],by="Sample") %>% 
  filter(egfr_status=="egfr_mut")
res.cut <- surv_cutpoint(dat2, time = "time", event = "status",
                         variables = "SRGN" )
summary(res.cut)
res.cat <- surv_categorize(res.cut)
dat2$SRGN_group<-res.cat$SRGN
dat<-rbind(dat1,dat2)
dat$group<-ifelse(dat$egfr_status=="egfr_wt",ifelse(dat$SRGN_group=="low","group1","group2"),ifelse(dat$SRGN_group=="low","group3","group4"))
load('02Cluster_IM.Rda')
IMM<-data.table::fread("./02Immune_Feature.txt")
dat_im<-dat %>% inner_join(IMM,by="Sample")
exp<-exp %>% as.data.frame()
PDL1<-data.frame(PDCD1=exp["PDCD1",dat_im$Sample] )
IM<-IM[dat_im$Sample,]
dat_im<-cbind(dat_im,IM,PDL1)

table(dat$group)
group<-factor(dat$group)
design <- model.matrix(~0+group)
colnames(design)<-substring(colnames(design),6,)
rownames(design) <- dat$Sample
es_wt<-es_all[,dat$Sample]
contrast.matrix<-makeContrasts("group1-group2","group1-group3","group3-group4","group2-group4",
                               levels = design)
vfit1 <- lmFit(es_wt,design)
vfit1 <- contrasts.fit(vfit1, contrast.matrix) 
efit1<- eBayes(vfit1)
res <- decideTests(efit1, p.value=0.05,lfc=0.5)
summary(res)
DEGs1<-topTable(efit1,coef=1,n=Inf,p.value=0.05,adjust="BH",lfc=0.5)%>%arrange(desc(abs(logFC)))
DEGs2<-topTable(efit1,coef=2,n=Inf,p.value=0.05,adjust="BH",lfc=0.5)%>%arrange(desc(abs(logFC)))
DEGs3<-topTable(efit1,coef=3,n=Inf,p.value=0.05,adjust="BH",lfc=0.5)%>%arrange(desc(abs(logFC)))
DEGs4<-topTable(efit1,coef=4,n=Inf,p.value=0.05,adjust="BH",lfc=0.5)%>%arrange(desc(abs(logFC)))
cat<-unique(c(rownames(DEGs1)[1:3],rownames(DEGs2)[1:3],rownames(DEGs3)[1:3],rownames(DEGs4)[1:3]))
cat<-cat[!(is.na(cat))]
dat<-dat %>% arrange(group)
res<-es_all[cat,dat$Sample]

group<-dat[,c('group','Sample')]%>%
  left_join(target[,c("Sample","Cluster")],by="Sample") %>% column_to_rownames("Sample")
table(group$group)

anno_colors<-list(Cluster=c(Cluster1 ="#48466d",Cluster2="#3d84a8",Cluster3="#46cdcf"),
                  group=brewer.pal(4,"Set2"))
names(anno_colors$group)<-c("group1","group2","group3","group4")
library(pheatmap)
library(RColorBrewer)
rownames(res)<-tolower(rownames(res))
rownames(res)<-substring(rownames(res),6)
pdf("./EGFR_SRGN_heatmap.pdf")
pheatmap(res, annotation_col = group,annotation_colors = anno_colors,
         color = pal_gsea("default", n = 100, reverse = F)(100),
         cluster_rows  =F,
         cluster_cols = F,
         fontsize=10,
         fontsize_row=12,
         cellwidth =0.5,
         cellheight =18,
         show_colnames=F,
         show_rownames=T,
         fontsize_col=10)
dev.off()
### Figure 8 A-B Immunoinhibitor & stimulator----
dir.create("./EGFR_SRGN/")
my_comparisons<-list(c('group1','group2'),c('group3','group4'),c('group1','group3'),c('group2','group4'))
color<-brewer.pal(4,"Set2")
names(color)<-c("group1","group2","group3","group4")
colnames(dat_im)
pdf("./EGFR_SRGN/Immunoinhibitor.pdf",height=4,width=4)
ggplot(dat_im ,aes(group,-Immunoinhibitor,fill=group))+
  geom_boxplot()+
  theme_classic()+
  labs(x = "", y = "Immunoinhibitor")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=color)+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,label="p.signif",method="wilcox.test")+  
  stat_compare_means(method="kruskal.test")
dev.off()



## EGFR/SRGN expression Figure 7 L----
load("./02.TPM&Targets.Rda")
Targets<-rbind(dat1,dat2)
sample<-Targets %>% 
  #filter(exp_group=="low") %>% 
  # filter(egfr_status=="egfr_wt") %>% 
  select("Sample")
egfr_srgn<-log2(TPM_sym[c("EGFR","SRGN"),sample$Sample] +1)%>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(Targets[,c("Cluster","egfr_status","methy_group","exp_group","Sample")],by="Sample")

cg02851793<-beta_norm["cg02851793",sample$Sample]
t.test(SRGN~egfr_status,data=egfr_srgn)
high
Welch Two Sample t-test
data:  SRGN by egfr_status
t = -7.0334, df = 62.21, p-value = 1.853e-09
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -1.3042109 -0.7269664
sample estimates:
  mean in group egfr_mut  mean in group egfr_wt 
8.593541               9.609129 

low
Welch Two Sample t-test

data:  SRGN by egfr_status
t = -4.9718, df = 32.795, p-value = 2.034e-05
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.9687990 -0.4060607
sample estimates:
  mean in group egfr_mut  mean in group egfr_wt 
6.926177               7.613607 






save(dat1,dat2,file="EGFR_SRGN_Stratified")


Targets$group<-paste0(Targets$egfr_status,"_",Targets$exp_group)
table(Targets$group)
levels(Targets$group)
Targets$group[Targets$group=="egfr_wt_low"]<-"group1"
Targets$group[Targets$group=="egfr_wt_high"]<-"group2"
Targets$group[Targets$group=="egfr_mut_low"] <-"group3"
Targets$group[Targets$group=="egfr_mut_high"]<-"group4"


table(Targets$group)
library(RColorBrewer)
color<-brewer.pal(4,"Set2")
names(color)<-c("group1","group2","group3","group4")
colnames(dat_im)
my_comparisons<-list(c("group1,group3"),c("group2","group4"))
pdf("./EGFR_SRGN/SRGN_revised_byegfr.pdf",height=4,width=5)
ggplot(Targets ,aes(egfr_status,SRGN,fill=group))+
  geom_violin()+
  theme_classic()+
  facet_wrap("exp_group",ncol=2)+
  labs(x = "", y = "Expression of SRGN")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=color)+
  scale_y_continuous()+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12),
        strip.background=element_blank(),
        strip.text = element_text(colour = "black",size=12))
#stat_compare_means(label="p.signif",method="t.test")+  
#stat_compare_means()
dev.off()

### Correlation between EGFR and SRGN Figure7 M-----
cor.test(x=dat[dat$egfr_status=="egfr_mut",'SRGN'],y=dat[dat$egfr_status=="egfr_mut",'CYT'],method="spearman")
cor.test(x=dat[dat$egfr_status=="egfr_wt",'SRGN'],y=dat[dat$egfr_status=="egfr_wt",'CYT'],method="spearman")
Targets<-Targets %>% 
  inner_join(egfr_srgn[,c("Sample","EGFR")],by="Sample") 
pdf("./EGFR_SRGN/CoreelationofSRGNEGFR_revised.pdf",height=5,width=4)
ggplot(Targets,aes(x=EGFR,y=SRGN))+
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

cor.test(Targets$EGFR[Targets$exp_group=="high"],Targets$SRGN[Targets$exp_group=="high"],method="spearman")

cor.test(Targets$EGFR[Targets$exp_group=="low"],Targets$SRGN[Targets$exp_group=="low"],method="spearman")

high
cor:  rho -0.4638282 p-value = 1.44e-05
low 
cor:  rho 0.04860008  p-value = 0.3597



