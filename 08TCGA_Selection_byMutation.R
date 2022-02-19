## Mutation information-----
rm(list=ls())
load("02Cluster_IM.Rda")
load("./02.TPM&Targets.Rda")
exp0<-log2(TPM_sym+1)
rm(Targets,TPM,TPM_sym)
target<-data.table::fread("./Cluster.txt")%>%arrange(Cluster)
CYT<-exp0[c("GZMA","PRF1"),]%>%
  t()%>%as.data.frame()%>%
  mutate(CYT=1/2*(log(PRF1)+log(GZMA)))
IM<-t(res_gsva)[target$Sample,]
IM<-cbind(IM[target$Sample,],CYT=CYT[target$Sample,"CYT"])
rm(res_gsva,IPS,Type,CYT)
target$age_group[target$Age<65]<-"<65"
target$age_group[target$Age>=65]<-">=65"
table(target$Stage)
target$stage_group[target$Stage%in%c("Stage I","Stage II")]<-"StageI_II"
target$stage_group[target$Stage%in%c("Stage III","Stage IV")]<-"StageIII_IV"
snp<-data.table::fread("./CNV_SNP/luadsnp.txt")
Targets<-data.table::fread("./Cluster.txt")
colnames(Targets)[7]<-"Tumor_Sample_Barcode"
egfr<-snp[snp$Hugo_Symbol=="EGFR",c(1,2,5,6,7,9,10,16,35:39,40:43)] %>% 
  inner_join(Targets[,c(7,21)],by="Tumor_Sample_Barcode")
egfr_mut<-unique(egfr$Tumor_Sample_Barcode)
egfr_wt<-setdiff(target$Sample,egfr_mut)
target<-target %>% 
  inner_join(data.frame(Sample=c(egfr_mut,egfr_wt),
                        egfr_status=c(rep("egfr_mut",length(egfr_mut)),rep("egfr_wt",length(egfr_wt)))),
             by="Sample")

tp53<-snp[snp$Hugo_Symbol=="TP53",c(1,2,5,6,7,9,10,16,35:39,40:43)]%>%
  inner_join(Targets[,c(7,21)],by="Tumor_Sample_Barcode")
tp53_mut<-unique(tp53$Tumor_Sample_Barcode)
tp53_wt<-setdiff(target$Sample,tp53_mut)
target<-target %>% 
  inner_join(data.frame(Sample=c(tp53_mut,tp53_wt),
                        tp53_status=c(rep("tp53_mut",length(tp53_mut)),rep("tp53_wt",length(tp53_wt)))),
             by="Sample")
target$N_Stage[target$N_Stage%in%c('N1','N2','N3')]<-"N"
rm(snp,egfr_mut,egfr_wt,egfr,tp53,tp53_mut,tp53_wt,Targets)
write.table(target,"./Cluster_mut.txt",sep="\t",quote=F,col.names = T,row.names = F)


##survival group-----
library(survival)
rm(list=ls())
load("./02.TPM&Targets.Rda")
target<-data.table::fread("./Cluster_mut.txt")%>%arrange(Cluster) 
# egfr_mut  egfr_wt 
# 59      391 
exp0<-log2(TPM_sym+1)
rm(TPM,TPM_sym)
load("./04CHAMPline_import.Rda")
beta0<-beta_norm
rm(beta_norm,Myload)

load("/Users/huzixin/Data_Project/Methylation/00mIMg.rda")
mIMg<-result
rownames(exp0)<-str_replace_all(rownames(exp0),"-","_")
probe<-mIMg$Probe
gene<-unique(mIMg$gene)
gene<-str_replace_all(gene,"-","_")
##common codes------
table(target$egfr_status)
Clinic<-target%>%
  filter(!is.na(time)&!(time==0))%>%
  filter(KRAS_status=="KRAS_wt")

colnames(Clinic)[c(9,11)]<-c("time","status")
Surv_ob<-Surv(time=Clinic$time,event=Clinic$status) 
exp<-exp0[,Clinic$Sample]%>%
  dplyr::filter(rownames(exp0)%in%gene)%>%
  t()%>%as.data.frame()
beta<-beta0[,Clinic$Sample]%>%as.data.frame()%>%
  dplyr::filter(rownames(beta0)%in%probe)%>%
  t()%>%as.data.frame()


##Univariate Cox------
colnames(Clinic)
data<-cbind(t(exp0[gene,Clinic$Sample]),
            Clinic[,c(9,11,10)])
Unicox_result<-data.frame(matrix(NA,ncol=4,nrow=length(gene)))
for(i in 1:length(gene)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', gene[i],collapse="+")) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") 
  Unicox_result[i,]<-c(gene[i],HR,CI,PValue)}
colnames(Unicox_result)<-c("gene","HR","95%CI","P.val")
Unicox_result$HR<-as.numeric(Unicox_result$HR) 
Unicox_result$P.val<-as.numeric(Unicox_result$P.val)

data<-cbind(t(beta0[probe,Clinic$Sample]),
            Clinic[,c(9,11)])
unicox_result<-data.frame(matrix(NA,ncol=4,nrow=length(probe)))
for(i in 1:length(probe)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = probe[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', probe[i],collapse="+"))
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2],2)
  PValue <- round(GSum$coefficients[,5],3)
  CI <- paste0(round(GSum$conf.int[,3:4],2),collapse = "-") 
  unicox_result[i,]<-c(probe[i],HR,CI,PValue)}

colnames(unicox_result)<-c("Probe","HR","95%CI","P.val")
unicox_result$HR<-as.numeric(unicox_result$HR)
unicox_result$P.val<-as.numeric(unicox_result$P.val)
UniCox_result_tcga<-Unicox_result
UniCox_result_tcga_probe<-unicox_result
UniCox_result_tcga_probe<-UniCox_result_tcga_probe %>% 
  inner_join(imps_Pheno[,c("gene","Probe")],by="Probe")
unique(intersect(UniCox_result_tcga_probe$gene[UniCox_result_tcga_probe$P.val<0.05],
                 UniCox_result_tcga$gene[UniCox_result_tcga$P.val<0.05]))
save(UniCox_result_tcga_probe,UniCox_result_tcga,file="./tcga_Survival_Results/tcga_unicox2(egfr_wt).rda")
### Multicox--------
data<-cbind(t(exp0[gene,Clinic$Sample]),
            Clinic[,c(9,11,10,22,23)])

colnames(data)[89:91]
Surv_ob<-Surv(data$time, data$status)
result<-data.frame(matrix(NA,ncol=5,nrow=length(gene)))
for(i in 1:length(gene)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = gene[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(c(gene[i],colnames(data)[89:91]),collapse="+"))) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2][1],2)
  PValue <- round(GSum$coefficients[,5][1],3)
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result[i,]<-c(gene[i],HR,names(HR),CI,PValue)}
colnames(result)<-c("gene","HR","group","95%CI","P.val")
result$HR<-as.numeric(result$HR)
result$P.val<-as.numeric(result$P.val)

data<-cbind(t(beta0[probe,Clinic$Sample]),
            Clinic[,c(9,11,10,22,23)])
colnames(data)[171:173]
Surv_ob<-Surv(data$time, data$status)
result2<-data.frame(matrix(NA,ncol=5,nrow=length(probe)))
for(i in 1:length(probe)){
  res.cut <- surv_cutpoint(data, time = "time", event = "status",variables = probe[i])
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  data[,i]<-res.cat[,3]
  FML <- as.formula(paste0('Surv_ob~', paste0(c(probe[i],colnames(data)[20:22]),collapse="+"))) 
  GCox <- coxph(FML,data=data) 
  GSum <- summary(GCox) 
  HR <- round(GSum$coefficients[,2][1],2)
  PValue <- round(GSum$coefficients[,5][1],3)
  CI <- paste0(round(GSum$conf.int[,3:4][1,],2),collapse = "-") 
  result2[i,]<-c(probe[i],HR,names(HR),CI,PValue)}
colnames(result2)<-c("Probe","HR","group","95%CI","P.val")
result2$HR<-as.numeric(result2$HR)
result2$P.val<-as.numeric(result2$P.val)
result_tcga<-result
result_tcga_probe<-result2
save(result_tcga,result_tcga_probe,
     file="./tcga_Survival_Results/tcga_multicox4(kras_wt).rda")



egfr_wt_tcga<-imps_Pheno %>% 
  filter(Probe%in%result2$Probe[result2$P.val<0.05]) %>% 
  filter(gene%in% result$gene[result$P.val<0.05]) %>% 
  dplyr::select(c('Probe','gene'))
tp53_wt_tcga<-imps_Pheno %>% 
  filter(Probe%in%result2$Probe[result2$P.val<0.05]) %>% 
  filter(gene%in% result$gene[result$P.val<0.05]) %>% 
  dplyr::select(c('Probe','gene'))
kras_wt_tcga<-imps_Pheno %>% 
  filter(Probe%in%result2$Probe[result2$P.val<0.05]) %>% 
  filter(gene%in% result$gene[result$P.val<0.05]) %>% 
  dplyr::select(c('Probe','gene'))

