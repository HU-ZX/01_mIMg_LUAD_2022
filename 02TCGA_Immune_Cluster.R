## Figure 2  
#### 04 Immunophenogram and Cluster----
### 0401 Cluster----
rm(list=ls())
load("02.TPM&Targets.Rda")
library(GSVA)
library(glmnet) 
library(tidyverse)
exp<-TPM_sym%>%as.data.frame()
unique(TIL$Cell.type)

TIL<-read.table("/Volumes/T7/immune_tmb/immune_related_biomarkers/ssGSEA/Charoentong.txt",sep="\t",header=T)
TIL<-TIL[,1:2]
colnames(TIL)<-c("GENE","CATEGORY")
IMM_GENES <- data.table::fread("/Volumes/T7/Methylation/GEO/IMM_GENES.txt")
unique(IMM_GENES$CATEGORY)
CM<-rbind(TIL,IMM_GENES)


cellMarker<- lapply(split(CM,CM$CATEGORY), function(x){
  dd = x$GENE
  unique(dd)
})

res_gsva <- gsva(as.matrix(exp),cellMarker, method = "ssgsea",kcdf="Gaussian")

IM<-res_gsva%>%t()%>%as.data.frame()%>%
  dplyr::select(c(2,3,9,10,4,19,28,21,22,24,25,30,32,15,16))
colnames(IM)<-c("aCD4","aCD8","emCD4","emCD8","aDC",
                "MDSC","Treg","MHCI","MHCII","NK","NKT","Th1","Th2","Immunoinhibitor","Immunostimulator")
IM<-IM%>%
  mutate(Effector=aCD4+aCD8+emCD4+emCD8+aDC+NK+NKT,
         Suppressor=-(MDSC+Treg),
         MHC=MHCI+MHCII,
         Immunoinhibitor=-Immunoinhibitor,
         Th2=-Th2)%>%
  dplyr::select(Effector,Suppressor,MHC,Th1,Th2,Immunoinhibitor,Immunostimulator)
g_im<-data.frame(as.matrix(IM))%>%dist()%>%hclust()

IM<-as.matrix(IM)
d<- sweep(t(IM),1, apply(IM,1,median,na.rm=T))
library(ConsensusClusterPlus)
res<-ConsensusClusterPlus(d=t(d),maxK=6,reps=1000,pItem=0.8,
                          pFeature=0.8,
                          clusterAlg = "km",distance="euclidean")
icl <- calcICL(res,title="test_Clusterplus",plot="png")

y=cutree(g_im,3)   
library(sparcl)
ColorDendrogram(g_im, y = y, labels = names(y), branchlength = 0.3,xlab=" ",sub=" ",main = " ")
Type<-as.data.frame(y)%>%rownames_to_column("Sample")%>%
  mutate(Cluster=paste0("Cluster",y))%>%
  dplyr::select(-2)%>%
  arrange(Cluster)
table(Type$Cluster)
# Cluster1 Cluster2 Cluster3 
# 247      166       37 

Targets<-Targets%>%inner_join(Type,by="Sample")
write.table(Targets,file="./Cluster.txt",col.names = T,row.names = F,quote=F,sep="\t")
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
library(ggplot2)
library(ggpubr)
dir.create("./Cluster_Phenotype_Plot/")
###  CYT Figure 2 B----
load("02.TPM&Targets.Rda")
exp<-TPM_sym%>%as.data.frame()
Targets<-data.table::fread("./Cluster.txt")
CYT<-exp[c("GZMA","PRF1"),]%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("Sample")
CYT<-CYT%>%
  mutate(CYT=1/2*(log(CYT$PRF1)+log(CYT$GZMA)))
CYT<-CYT%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")

with(CYT,tapply(CYT,Cluster,shapiro.test))
bartlett.test(CYT~Cluster,data = CYT) 
with(CYT,tapply(CYT,Cluster,mean))
#Cluster1 Cluster2 Cluster3 
#2.661860 3.730718 1.848335
pdf("./Cluster_Phenotype_Plot/01CYT.pdf",width=4,height=4)
ggplot(CYT,aes(Cluster,CYT,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  labs(x = "", y = "Cytolytic activity (CYT)")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'))+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,method = "t.test",label="p.signif")+  
  stat_compare_means(method="anova",label.x=0.6,label.y=7)
dev.off()

## IM_score Figure 2 C-----
IM_score<-IM%>%
  as.data.frame() %>% 
  mutate(Im_score=Effector+Suppressor+MHC+Th1+Th2+Immunoinhibitor+Immunostimulator)%>%
  rownames_to_column("Sample")%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")
with(IM_score,tapply(Im_score,Cluster,shapiro.test))
bartlett.test(Im_score~Cluster,data = IM_score) 
pdf("./Cluster_Phenotype_Plot/TCGA_IMpS.pdf",width=4,height=4)
ggplot(IM_score,aes(Cluster,Im_score,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  labs(x = "", y = "Immunophenotype Score (IMpS)")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'))+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,label="p.signif")+  
  stat_compare_means(method="kruskal.test",label.x=0.6,label.y=7.8)
dev.off()

##PD1----
exp<-log2(exp+1)
PD1<-exp["CD274",]%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("Sample")%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")
with(PD1,tapply(CD274,Cluster,shapiro.test))
pdf("./Cluster_Phenotype_Plot/01PD1.pdf",width=4,height=4)
ggplot(PD1,aes(Cluster,CD274,fill=Cluster))+
  geom_violin()+
  theme_classic()+
  labs(x = "", y = "PD-L1 expression")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'))+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,label="p.signif",method="wilcox.test")+  
  stat_compare_means(method="kruskal.test",label.x=0.8,label.y=9.5)
dev.off()


### Heatmap Figure 2 A------
dat<-IM[,c("MHC","Th1","Effector","Immunostimulator",
           "Th2","Suppressor","Immunoinhibitor")]
Cluster<-Type%>%arrange(Cluster)
dat<-dat[Cluster$Sample,]%>%scale()

Cluster<-Cluster%>%column_to_rownames("Sample")
table(Cluster$Cluster)
anno_colors<-list(Cluster=c(Cluster1 ="#48466d",Cluster2="#3d84a8",Cluster3="#46cdcf"))
library(pheatmap)
pdf("./Cluster_Phenotype_Plot/heatmap.pdf",width=4,height=4)
pheatmap(dat, annotation_row = Cluster, annotation_colors = anno_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows  =F,
         cluster_cols = F,
         fontsize=10,
         fontsize_row=8,
         cellwidth =25,
         cellheight =0.4,
         show_colnames=T,
         show_rownames=F,
         fontsize_col=10)

dev.off()
quantile(dat)
save(res_gsva,CYT,Type,IM,IPS,file="02Cluster_IM.Rda")


### Beanplot Figure 2 D ----
rm(list=ls())
library(EpiDISH)
library(tidyverse)
library(ggpubr)
data(centEpiFibIC.m)
data(DummyBeta.m)
data(LiuDataSub.m)
data(centBloodSub.m)
load("./02.TPM&Targets.Rda")
load("./04CHAMPline_import.Rda")
Targets<-data.table::fread("Cluster_mut.txt")
frac.m <- hepidish(beta.m = beta_norm, ref1.m = centEpiFibIC.m, ref2.m = centDHSbloodDMC.m, h.CT.idx = 3, method = 'RPC')
frac.m <-frac.m %>%as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")
frac.m$CD3<-frac.m$CD4T+frac.m$CD8T
with(frac.m,tapply(CD3,Cluster,shapiro.test))
with(frac.m,tapply(Epi,Cluster,shapiro.test))
bartlett.test(Epi~Cluster,data =frac.m) 
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
dat<-frac.m %>% 
  arrange(Cluster)%>% 
  column_to_rownames("Sample")%>% 
  select(-"Cluster") %>% 
  select("CD3","Fib","Epi")
dat_c<-frac.m %>% 
  arrange(Cluster)%>% 
  column_to_rownames("Sample")%>% 
  select("CD3","Epi","Cluster") %>% 
  reshape2::melt(id=c("Cluster"))
dat_c$Cells<-as.character(dat_c$Cells)
colnames(dat_c)[2:3]<-c("Cells","Score")


kruskal.test(CD3~Cluster,frac.m)
# data:  CD3 by Cluster
# Kruskal-Wallis chi-squared = 160.76,
# df = 2, p-value < 2.2e-16
#两两比较
with(frac.m,pairwise.wilcox.test(CD3,Cluster,exact=F))
#          Cluster1 Cluster2
# Cluster2 < 2e-16  -       
# Cluster3 2.3e-07  < 2e-16 
kruskal.test(Epi~Cluster,frac.m)
# Kruskal-Wallis chi-squared = 142.69,
# df = 2, p-value < 2.2e-16
with(frac.m,pairwise.wilcox.test(Epi,Cluster,exact=F))
#         Cluster1 Cluster2
# Cluster2 < 2e-16  -       
# Cluster3 0.00076  8.3e-14 

pdf("./Cluster_Phenotype_Plot/03CD3_Epi_byEpidish.pdf",width=5,height=5)
beanplot(Score ~ Cells*Cluster,dat_c,
         side = "b",
         col = list("#4DBBD5FF","#00A087FF"), 
         border = c("#2C738EFF","#4FC46AFF"), 
         main = "p-value < 2.2e-16 ",
         xlab = "", ylab = "Fraction of cells",ylim = c(0,1))
legend("topleft",bty="n",c("CD3+T cells", "Epithelial cells"),
       fill = c("#4DBBD5FF","#00A087FF"))
dev.off()


ggplot(frac.m,aes(Cluster,CD3,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  labs(x = "", y = "Cytolytic activity (CYT)")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position="none")+
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'),labels=c('Cluster 1','Cluster 2','Cluster 3'))+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,method="t.test",label="p.signif")+  
  stat_compare_means()


