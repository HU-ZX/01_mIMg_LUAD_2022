### Data Import-----
library(GEOquery)
library(tidyverse)
library(ggpubr)
rm(list=ls())
files<-list.files(".", full.names = T)
methy_set<-getGEO(GEO="GSE119144",file=files[4],getGPL=F)
exp_set<-getGEO(GEO="GSE135222",file=files[6],getGPL=F)
target<-pData(exp_set)
# "Supplementary_files_format_and_content: gene expression matrix (tab-delimited text file; TPM unit)"
colnames(target)
target<-target[,c(1,41:44)]
target$title<-substring(target$title,7,)
colnames(target)<-c("ID","Age","Gender","PFS_time","PFS")
target2<-pData(methy_set)
colnames(target2)
target2<-target2[,1:2]
colnames(target2)<-c("ID","geo_accession")
target2$ID<-str_extract(target2$ID,"\\d+")
Table_S6 <- readxl::read_excel("~/Data_Project/Methylation/ICIs_Response/GSE126043&GSE126044/Supplementary_Table_S1-6/Supplementary Table S6.xlsx")
colnames(Table_S6)<-Table_S6[1,]
Table_S6<-Table_S6[-1,]
colnames(Table_S6)[c(1,4)]<-c("geo_accession","ID")
target2<-target2 %>% 
  inner_join(Table_S6[,-1],by="ID")
methy<-exprs(methy_set)[,target2$geo_accession]
colnames(methy)<-target2$ID
count<-data.table::fread(files[5])
load("/Volumes/T7/Index_Gencode_v35/anno.Rdata")

count<-count%>%
  mutate(gene_id=substring(count$gene_id,1,15))%>%
  column_to_rownames("gene_id")
colnames(count)<-str_extract(colnames(count),"\\d+")
count<-count[,target$ID]
target<-target %>% 
  inner_join(Table_S6[,-1],by="ID")


TPM<-count

identical(colnames(TPM),target$ID)
TPM_sym<-TPM%>%
  mutate(ENSG=rownames(TPM))%>%
  inner_join(anno[,4:5],by="ENSG")%>%
  distinct(Symbol,.keep_all = T)%>%
  column_to_rownames("Symbol")%>%
  dplyr::select(-"ENSG")
exp<-TPM_sym
save(methy,exp,target,target2,file="./ICIs_Response2/data2.rda")
### Figure 3 D all_methylation-----
rm(list=ls())
load("./ICIs_Response2/data2.rda")
all_methy<-apply(methy,2,mean)
all_methy<-data.frame(ID=names(all_methy),methy=all_methy) 
all_methy<-all_methy%>%
  inner_join(target2[,c("ID","Responsiveness")],by="ID") %>% 
  arrange("Responsiveness")%>% 
  mutate(Responsiveness=ifelse(Responsiveness=="DCB","Responders","Non-responders"))
with(all_methy,tapply(methy,Responsiveness,quantile))
wilcox.test(methy~Responsiveness,data=all_methy,paired=F,conf.int=T,correct=F,conf.level=0.95)

pdf("./ICIs_Response2/data2_Methylation_all.pdf",width=4,height=4)
ggplot(all_methy,aes(Responsiveness,methy))+
  geom_boxplot()+
  theme_bw()+
  labs(x = "", y = "Average methylation level")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank(),
        legend.position="none") +
  scale_y_continuous()+ 
  stat_compare_means(method="wilcox.test")

dev.off()
### Figure5 F-I Validation -----
rm(list=ls())
library(survival)
library(survminer) 
load("./ICIs_Response2/data2.rda")
target2<-target2 %>% 
  mutate(Responsiveness=ifelse(Responsiveness=="DCB","Responders","Non-responders"))
target<-target %>% 
  mutate(Responsiveness=ifelse(Responsiveness=="DCB","Responders","Non-responders"))
colnames(target)

load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%probeset1$Probe,"gene"])

gene<-unique(probe_exp_meth_cr$gene)
probe<-unique(probe_exp_meth_cr$Probe)
gene<-unique(geneset1)
probe<-intersect(probeset1$Probe,rownames(methy))
gene<-unique(imps_Pheno$gene)
probe<-unique(imps_Pheno$Probe)
methy<-methy[,target2$ID]
diff<-data.frame()
j=1
probe<-intersect(rownames(methy),probe)
dir.create("./ICIs_Response2/Response_Plot/")
for(i in 1:160){
  dat<-data.frame(probe=methy[probe[i],]) %>% 
    rownames_to_column("ID") %>% 
    inner_join(target2,by="ID")
  test<-wilcox.test(probe~Responsiveness,data=dat,paired=F,conf.int=T,correct=F,conf.level=0.95)
  m1<-mean(dat$probe[dat$Responsiveness=="Non-responders"])
  m2<-mean(dat$probe[dat$Responsiveness=="Responders"])
  probe_beta<-imps_Pheno$tcga_imps[imps_Pheno$Probe==probe[i]]
  meidan<-tapply(dat$probe,dat$Responsiveness,quantile)
  meidan_r<-paste0(round(meidan[[2]][3],2)," (",round(meidan[[2]][2],2),",",round(meidan[[2]][4],2),")")
  meidan_n<-paste0(round(meidan[[1]][3],2)," (",round(meidan[[1]][2],2),",",round(meidan[[1]][4],2),")")
  conf<-round(test$conf.int,2)
  d<-paste0(round(test$estimate,2)," (",conf[1],"-",conf[2],")")
  p.val<-round(test$p.value,3)
  # pdf(paste0("./ICIs_Response2/Response_Plot/", 
  #            imps_Pheno$gene[imps_Pheno$Probe==probe[i]],
  #            "_",
  #            probe[i],".pdf"),height=4,width=4)
  # boxplot(probe~Responsiveness,data=dat,ylab="beta value",
  #         xlab="",
  #         main="")
  # mtext(paste0(probe[i],"/",
  #              imps_Pheno$gene[imps_Pheno$Probe==probe[i]],
  #              "\n",
  #              "Mann-Whitney U-test","\n","p.val = ",round(p.val,3)),
  #       cex=0.8,padj=-0.2,adj=0.1,font=3
  #       )
  # dev.off()
  diff[j,1]<-probe[i]
  diff[j,2]<-imps_Pheno$gene[imps_Pheno$Probe==probe[i]]
  diff[j,3]<-meidan_n
  diff[j,4]<-meidan_r
  diff[j,5]<-d
  diff[j,6]<-p.val
  j=j+1
  if(((m2-m1)>0&probe_beta>0)|((m2-m1)<0&probe_beta<0)){
    if(p.val<0.05){
      print(probe[i])
    }
  }
}

colnames(diff)<-c("Probe","gene","Median (P25, P75) (Non-responder)","Median (P25, P75) (Responder)","d (95%CI)(Non-responder minus Responder)","P.val")
diff$probe_gene<-paste0(diff$Probe,"/",diff$gene)
diff1<-diff %>% 
  inner_join(imps_Pheno[,c(1,11,10,12)],by="Probe") %>% 
  arrange(gene)
colnames(imps_Pheno)
diff %>% 
  filter(Probe%in%probe_exp_meth_cr$Probe)
write.table(diff1,"./ICIs_Response2/validate_diff.txt",
            sep="\t",col.names = T,row.names = F,
            quote=F)

diff2<-data.frame()
j=1
gene<-intersect(rownames(exp),gene)
exp<-log2(exp+1)
dir.create("./ICIs_Response2/Response_Plot_gene/")
for(i in 1:86){
  dat<-exp[gene[i],]%>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    inner_join(target,by="ID")
  colnames(dat)[2]<-"gene"
  test<-wilcox.test(gene~Responsiveness,data=dat,paired=F,conf.int=T,correct=F,conf.level=0.95)
  m1<-mean(dat$gene[dat$Responsiveness=="Non-responders"])
  m2<-mean(dat$gene[dat$Responsiveness=="Responders"])
  gene_imps<-imps_Pheno$tcga_imps_g[imps_Pheno$gene==gene[i]]
  meidan<-tapply(dat$gene,dat$Responsiveness,quantile)
  meidan_r<-paste0(round(meidan[[2]][3],2)," (",round(meidan[[2]][2],2),",",round(meidan[[2]][4],2),")")
  meidan_n<-paste0(round(meidan[[1]][3],2)," (",round(meidan[[1]][2],2),",",round(meidan[[1]][4],2),")")
  conf<-round(test$conf.int,2)
  d<-paste0(round(test$estimate,2)," (",conf[1],"-",conf[2],")")
  p.val<-round(test$p.value,3)
  pdf(paste0("./ICIs_Response2/Response_Plot_gene/", 
             gene[i],
             ".pdf"),height=4,width=4)
  boxplot(gene~Responsiveness,data=dat,ylab="beta value",
          xlab="",
          main="")
  mtext(paste0(gene[i],
               "\n",
               "Mann-Whitney U-test","\n","p.val = ",round(p.val,3)),
        cex=0.8,padj=-0.2,adj=0.1,font=3
  )
  dev.off()
  diff2[j,1]<-gene[i]
  diff2[j,2]<-meidan_n
  diff2[j,3]<-meidan_r
  diff2[j,4]<-d
  diff2[j,5]<-p.val
  j=j+1
  if(((m2-m1)>0&gene_imps>0)|((m2-m1)<0&gene_imps<0)){
    if(p.val<0.05){
      print(gene[i])
    }
  }
}
colnames(diff2)<-c("gene","Median (P25, P75) (Non-responder)","Median (P25, P75) (Responder)","d (95%CI)(Non-responder minus Responder)","P.val")
# [1] "RASSF4"
# [1] "FASLG"
# [1] "SRGN"
# [1] "SLAMF8"
# [1] "CXCR6"
colnames(imps_Pheno)
imps<-imps_Pheno[,c(2,14,13,15)] %>% distinct(gene,.keep_all = T)
diff3<-diff2 %>%arrange(gene) %>% 
  left_join(imps,by="gene") 
write.table(diff3,"./ICIs_Response2/validate_diff_gene.txt",
            sep="\t",col.names = T,row.names = F,
            quote=F)
# ACAP1  ARHGAP30 CD3D CD247 LCK LAX1 PSTPIP1 

vali_gene<-c('ACAP1',  'ARHGAP30', 'CD3D', 'CD247', 'LCK', 'LAX1', 'PSTPIP1' )
vali_probe<-probeset1$Probe[probeset1$gene %in%vali_gene]
j=1
diff<-data.frame()
for(i in 1:length(vali_probe)){
  dat<-data.frame(probe=methy[vali_probe[i],]) %>% 
    rownames_to_column("ID") %>% 
    inner_join(target2,by="ID")
  test<-kruskal.test(probe~Responsiveness,data=dat)
  p.val<-test$p.value
  if(p.val<0.05){
    diff[j,1]<-vali_probe[i]
    diff[j,2]<-probeset1$gene[probeset1$Probe==vali_probe[i]]
    diff[j,3]<-p.val
    j=j+1
  }
  boxplot(probe~Responsiveness,data=dat,ylab=vali_probe[i],main=paste0(probeset1$gene[probeset1$Probe==vali_probe[i]],"\n","p.val = ",p.val))
  dat<-data.frame(probe=methy[vali_probe[i],]) %>% 
    rownames_to_column("ID") %>% 
    inner_join(target,by="ID")
  dat$PFS<-as.numeric(dat$PFS)
  dat$PFS_time<-as.numeric(dat$PFS_time)
  res.cut <- surv_cutpoint(dat, time = "PFS_time", event = "PFS",variables = "probe")
  summary(res.cut)
  res.cat <- surv_categorize(res.cut)
  Surv_ob<-Surv(res.cat$PFS_time, res.cat$PFS)
  fit <- survfit(Surv_ob~probe, data = res.cat)
  model_dif<-survdiff(Surv_ob~probe, data = res.cat)
  pValue<- 1-pchisq(model_dif$chisq,df=1)
  plot(fit, 
       lwd=2,
       col=c("red","blue"),
       xlab="Time (month)",
       #mark.time=T,
       ylab="Survival rate",
       main=paste0(vali_probe[i],"(p.val =", pValue ,")",sep=" "))
  legend("topright", 
         c("High","Low"), 
         lwd=2, 
         col=c("red","blue")
  )
}

colnames(diff)<-c("probe","gene","p.val")
diff$probe
vali_probe2<-c("cg26227523" ,"cg07728874", "cg24841244")
dat<-data.frame(probe=methy["cg25671438",]) %>% 
  rownames_to_column("ID") %>% 
  inner_join(target2,by="ID")
test<-kruskal.test(probe~Responsiveness,data=dat)
p.val<-test$p.value
boxplot(probe~Responsiveness,data=dat,ylab=probe[i],main=paste0(probeset1$gene[probeset1$Probe==probe[i]],"\n","p.val = ",p.val))


dat<-t(methy[vali_probe,target2$ID]) %>% 
  as.data.frame() 
signature<-cbind(signature=apply(dat,1,sum),target2)
test<-kruskal.test(signature~Responsiveness,data=signature)
p.val<-test$p.value
boxplot(signature~Responsiveness,data=signature)
data<-target %>% 
  inner_join(signature[,c("signature","ID")],by="ID")

cor.test(as.numeric(exp["ACAP1",target$ID]),as.numeric(methy["cg25671438",target$ID]))
data<-data.frame(t(exp["ACAP1",target$ID]),
                 probe=methy["cg25671438",target$ID],
                 target)
colnames(data)[1]<-"gene"
data$PFS<-as.numeric(data$PFS)
data$PFS_time<-as.numeric(data$PFS_time)
res.cut <- surv_cutpoint(data, time = "PFS_time", event = "PFS",variables = "signature")
summary(res.cut)
res.cat <- surv_categorize(res.cut)

Surv_ob<-Surv(res.cat$PFS_time, res.cat$PFS)
fit <- survfit(Surv_ob~signature, data = res.cat)
model_dif<-survdiff(Surv_ob~signature, data = res.cat)
pValue<- 1-pchisq(model_dif$chisq,df=1)

tiff(paste0("./tcga_Survival_Results/plot/imps_gene/",gene[i],"_",".tiff"))
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (month)",
     #mark.time=T,
     ylab="Survival rate",
     main=paste0("ACAP1","(p.val =", pValue ,")",sep=" "))
legend("topright", 
       c("High","Low"), 
       lwd=2, 
       col=c("red","blue")
)
dev.off()
table(res.cat$probe)
res.cat2<-data.frame(probe=res.cat[,"probe"]) 
rownames(res.cat2)<-rownames(res.cat)
res.cat2<-res.cat2%>% 
  rownames_to_column("ID") %>% 
  inner_join(target,by="ID") %>% 
  mutate(new=paste0(probe,Responsiveness)) %>% 
  add_count(new) %>% 
  distinct(new,.keep_all = T) %>% 
  dplyr::select("Responsiveness","probe","n") %>% 
  pivot_wider(names_from = "probe",values_from = "n") %>% 
  column_to_rownames("Responsiveness")
fisher.test(res.cat2)

if(pValue<0.05){
  diff[j]<-gene[i]
  j=j+1
}

















### Survival




###
dat1<-exp[c("PDCD1","LCK","CD247","PSTPIP1"),target$ID] %>% t() %>% as.data.frame()
cor.test(dat1$PDCD1,dat1$CD247)
dat2<-methy[c('cg09032544', 'cg07786657', 'cg11683242', 'cg26227523'),target$ID] %>% t() %>% as.data.frame()
cor.test(dat1$PDCD1,dat2$cg26227523)
dat1<-dat1 %>% rownames_to_column("ID") %>% 
  inner_join(target,by="ID")
wilcox.test(CD247~Responsiveness,data=dat1)
boxplot(LCK~Responsiveness,data=dat1)





