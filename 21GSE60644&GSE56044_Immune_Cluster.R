###  Cluster----
library(sparcl)
load("./01Importdata.rda")
table(pda_co$histology)
norm<-norm[pda_co$title[pda_co$histology=="AC"]]
pda_ac<-pda_co[pda_co$histology=="AC",]

library(GSVA)
library(tidyverse)
exp<-norm%>%as.data.frame()


TIL<-read.table("/Volumes/T7/immune_tmb/immune_related_biomarkers/ssGSEA/Charoentong.txt",sep="\t",header=T)
TIL<-TIL[,1:2]
colnames(TIL)<-c("GENE","CATEGORY")
IMM_GENES <- data.table::fread("/Volumes/T7/Methylation/GEO/IMM_GENES.txt")
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

y=cutree(g_im,3)   
library(sparcl)
ColorDendrogram(g_im, y = y, labels = names(y), branchlength = 0.3,xlab=" ",sub=" ",main = " ")
Type<-as.data.frame(y)%>%rownames_to_column("title")%>%
  mutate(Cluster=paste0("Cluster",y))%>%
  dplyr::select(-2)%>%
  arrange(Cluster)
table(Type$Cluster)
# Cluster1 Cluster2 Cluster3 
#  36       31       11 
pda_ac<-pda_ac%>%inner_join(Type,by="title")%>%rename("Sample"="title")
write.table(pda_ac,file="./Cluster.txt",col.names = T,row.names = F,quote=F,sep="\t")

dir.create("./GSE60644_GSE60645_Data/Cluster_Phenotype_Plot/")
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
library(ggplot2)
library(ggpubr)
?shapiro.test()

###  CYT----
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
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'),labels=c('Cluster 1','Cluster 2','Cluster 3'))+
  scale_y_continuous()+
  stat_compare_means(comparisons=my_comparisons,method="t.test",label="p.signif")+  
  stat_compare_means(method="anova",label.y=2.6)
dev.off()
# IMPS------
IM_score<-IM%>%mutate(Im_score=Effector+Suppressor+MHC+Th1+Th2+Immunoinhibitor+Immunostimulator)%>%
  rownames_to_column("Sample")%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")
with(IM_score,tapply(Im_score,Cluster,shapiro.test))
bartlett.test(Im_score~Cluster,data = IM_score) 

pdf("./Cluster_Phenotype_Plot/GSM60644_IMpS.pdf",width=4,height=4)
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
  stat_compare_means(method="kruskal.test",label.x=0.9,label.y=3.3)
dev.off()


## PDL1-----
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
  stat_compare_means(method="kruskal.test",label.x=0.8)
dev.off()



#heatmap-----
dat<-IM[,c("MHC","Th1","Effector","Immunostimulator",
           "Th2","Suppressor","Immunoinhibitor")]
Cluster<-Type%>%arrange(Cluster)
dat<-dat[Cluster$title,]%>%scale()
Cluster<-Cluster%>%column_to_rownames("title")
table(Cluster$Cluster)
anno_colors<-list(Cluster=c(Cluster1 ="#48466d",Cluster2="#3d84a8",Cluster3="#46cdcf"))
library(pheatmap)
pdf("./Cluster_Phenotype_Plot/GSE60644_Cluster_heatmap.pdf",height=4,width=4)
pheatmap(dat, annotation_row = Cluster, annotation_colors = anno_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows  =F,
         cluster_cols  =F,
         fontsize=8,
         fontsize_row=8,
         cellwidth =20,
         cellheight =2,
         show_colnames=T,
         show_rownames=F,
         fontsize_col=8)
dev.off()
quantile(dat)

save(res_gsva,file="./GSE60644_GSE60645_Data/02Cluster.Rda")

# Figure 2 L Beanplot-----
rm(list=ls())
load("01Importdata.rda")
load("02Methylation.Rda")
Target<-data.table::fread("./Cluster.txt")
beta<-myNorm[,Target$Sample]
exp<-norm[,Target$Sample]
library(EpiDISH)
library(tidyverse)
library(ggpubr)
data(centEpiFibIC.m)
data(DummyBeta.m)
data(LiuDataSub.m)
data(centBloodSub.m)
frac.m <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centDHSbloodDMC.m, h.CT.idx = 3, method = 'RPC')
frac.m <-frac.m %>%as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  inner_join(Target[,c("Sample","Cluster")],by="Sample")
frac.m$CD3<-frac.m$CD4T+frac.m$CD8T
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
colnames(dat_c)[2:3]<-c("Cells","Score")
dat_c$Cells<-as.character(dat_c$Cells)
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
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


with(frac.m,tapply(CD3,Cluster,shapiro.test))
with(frac.m,tapply(Epi,Cluster,shapiro.test))
bartlett.test(Epi~Cluster,data =frac.m) 

kruskal.test(CD3~Cluster,frac.m)
# Kruskal-Wallis chi-squared = 47.474, 
# df = 2, p-value= 4.91e-11
with(frac.m,pairwise.wilcox.test(CD3,Cluster,exact=F))
#          Cluster1 Cluster2
# Cluster2 1.3e-08  -       
# Cluster3 0.0011   2.7e-06 
kruskal.test(Epi~Cluster,frac.m)
# Kruskal-Wallis chi-squared = 43.398, 
# df = 2, p-value= 3.769e-10
with(frac.m,pairwise.wilcox.test(Epi,Cluster,exact=F))
#         Cluster1 Cluster2
#Cluster2 7.3e-09  -       
#Cluster3 0.036    1.8e-05 
library(beanplot)
pdf("./Cluster_Phenotype_Plot/03CD3_Epi_byEpidish.pdf",width=5,height=5)
beanplot(Score ~ Cells*Cluster,dat_c,
         side = "b",
         col = list("#4DBBD5FF","#00A087FF"), 
         border = c("#2C738EFF","#4FC46AFF"), 
         main = "p-value = 4.91e-11",
         xlab = "", ylab = "Fraction of cells",ylim=c(0,1))
legend("topleft",bty="n",c("CD3+T cells", "Epithelial cells"),
       fill = c("#4DBBD5FF","#00A087FF"))
dev.off()
