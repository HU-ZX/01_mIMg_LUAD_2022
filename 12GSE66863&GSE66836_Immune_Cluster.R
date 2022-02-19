### Cluster----
rm(list=ls())
load("./01Importexpression.Rda")
library(GSVA)
library(tidyverse)
exp<-expr_norm%>%column_to_rownames("Gene")%>%as.data.frame()
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
IM<-t(res_gsva)
rm(IMM_GENES,cellMarker,CM,TIL)
save(IM,file="./02Immunecells_gse66863.rda")

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
Type<-as.data.frame(y)%>%rownames_to_column("Sample")%>%
  mutate(Cluster=paste0("Cluster",y))%>%
  dplyr::select(-2)%>%
  arrange(Cluster)
table(Type$Cluster)
# Cluster1 Cluster2 Cluster3 
#  77       21       23 
Type$Cluster[Type$Cluster%in%c("Cluster2","Cluster3")]<-ifelse(Type$Cluster[Type$Cluster%in%c("Cluster2","Cluster3")]=="Cluster2","Cluster3","Cluster2")
pda<-pda%>%mutate(Sample=geo_accession)%>%inner_join(Type,by="Sample")
Targets<-pda
write.table(pda,file="./Cluster.txt",col.names = T,row.names = F,quote=F,sep="\t")

my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
library(ggplot2)
library(ggpubr)
?shapiro.test()
### Figure 2 F CYT----
Targets<-data.table::fread("./Cluster.txt")%>%rename("Sample"="geo_accession_e")
CYT<-exp[c("GZMA","PRF1"),]%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("Sample")
CYT<-CYT%>%
  mutate(CYT=1/2*(log(CYT$PRF1)+log(CYT$GZMA)))
CYT<-CYT%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")

with(CYT,tapply(CYT,Cluster,shapiro.test))
bartlett.test(CYT~Cluster,data = CYT) 
dir.create("./Cluster_Phenotype_Plot/")
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
  stat_compare_means(method="anova",label.y=2.5)
dev.off()

## Figure 2 G IM_score----
IM_score<-IM%>%mutate(Im_score=Effector+Suppressor+MHC+Th1+Th2+Immunoinhibitor+Immunostimulator)%>%
  rownames_to_column("Sample")%>%
  inner_join(Targets[,c("Sample","Cluster")],by="Sample")
with(IM_score,tapply(Im_score,Cluster,shapiro.test))
bartlett.test(Im_score~Cluster,data = IM_score)
pdf("./Cluster_Phenotype_Plot/GSE66863_IMpS.pdf",width=4,height=4)
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
  stat_compare_means(method="kruskal.test",label.x=0.6)
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
  stat_compare_means(method="kruskal.test",label.x=0.8)dev.off()
dev.off()


# Figure 2 E Heatmap----
dat<-IM[,c("MHC","Th1","Effector","Immunostimulator",
           "Th2","Suppressor","Immunoinhibitor")]
Cluster<-Targets[,c("Sample","Cluster")]%>%arrange(Cluster)
dat<-dat[Cluster$Sample,]%>%scale()
Cluster<-Cluster%>%column_to_rownames("Sample")
table(Cluster$Cluster)
anno_colors<-list(Cluster=c(Cluster1 ="#48466d",Cluster2="#3d84a8",Cluster3="#46cdcf"))
library(pheatmap)
pdf("./Cluster_Phenotype_Plot/GSE66863_Cluster_heatmap.pdf",height=4,width=4)
pheatmap(dat, annotation_row = Cluster, annotation_colors = anno_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_rows  =F,
         cluster_cols  =F,
         fontsize=8,
         fontsize_row=8,
         cellwidth =25,
         cellheight =1.3,
         show_colnames=T,
         show_rownames=F,
         fontsize_col=8)
dev.off()

save(IM_score,file="./Cluster_Phenotype_Plot/gse66863Immunescore.rda")



## Figure 2 H Beanplot----
rm(list=ls())
library(EpiDISH)
data(centEpiFibIC.m)
data(DummyBeta.m)
data(LiuDataSub.m)
data(centBloodSub.m)
load("./01Importexpression.Rda")
load("./02Methylation_byChAMP.Rda")
Target<-data.table::fread("./Cluster.txt")
pD<-Target %>% 
  column_to_rownames("geo_accession_m")
pD<-pD[colnames(beta_norm),]

exp<-expr_norm[,pD$geo_accession_e]
rownames(exp)<-expr_norm$Gene
out.l <- epidish(beta.m = DummyBeta.m, ref.m = centEpiFibIC.m, method = "RPC") 
out.l$estF
dim(out.l$dataREF)
frac.m <- epidish(beta.m = beta_norm, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
boxplot(BloodFrac.m)
frac.m <- hepidish(beta.m = beta_norm, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
frac.m <-frac.m %>%as.data.frame() %>% 
  rownames_to_column("geo_accession_m") %>% 
  inner_join(Target[,c("geo_accession_m","Cluster")],by="geo_accession_m")

frac.m$CD3<-frac.m$CD4T+frac.m$CD8T
pdf("./Cluster_Phenotype_Plot/01CYT.pdf",width=4,height=4)
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
dev.off()







rm(list=ls())
load("./02Methylation_byChAMP.Rda")
Target<-data.table::fread("./Cluster.txt")
pD<-Target %>% 
  column_to_rownames("geo_accession_m")
pD<-pD[colnames(beta_norm),]
library(EpiDISH)
library(tidyverse)
library(ggpubr)
data(centEpiFibIC.m)
data(DummyBeta.m)
data(LiuDataSub.m)
data(centBloodSub.m)
frac.m <- hepidish(beta.m = beta_norm, ref1.m = centEpiFibIC.m, ref2.m = centDHSbloodDMC.m, h.CT.idx = 3, method = 'RPC')
frac.m <-frac.m %>%as.data.frame() %>% 
  rownames_to_column("geo_accession_m") %>% 
  inner_join(Target[,c("geo_accession_m","Cluster")],by="geo_accession_m")
frac.m$CD3<-frac.m$CD4T+frac.m$CD8T
dat<-frac.m %>% 
  arrange(Cluster)%>% 
  column_to_rownames("geo_accession_m")%>% 
  select(-"Cluster") %>% 
  select("CD3","Fib","Epi")
dat_c<-frac.m %>% 
  arrange(Cluster)%>% 
  column_to_rownames("geo_accession_m")%>% 
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
# Kruskal-Wallis chi-squared =29.471,
# df = 2, p-value =3.985e-07
with(frac.m,pairwise.wilcox.test(CD3,Cluster,exact=F))
#          Cluster1 Cluster2
# Cluster2 0.00016  -       
# Cluster3 0.00076  1.9e-05 
kruskal.test(Epi~Cluster,frac.m)
# Kruskal-Wallis chi-squared = 26.031, 
# df =2, p-value = 2.226e-06
with(frac.m,pairwise.wilcox.test(Epi,Cluster,exact=F))
#          Cluster1 Cluster2
# Cluster2 0.00192  -       
# Cluster3 0.00054  2.7e-05 
library(beanplot)
pdf("./Cluster_Phenotype_Plot/03CD3_Epi_byEpidish.pdf",width=5,height=5)
beanplot(Score ~ Cells*Cluster,dat_c,
         side = "b",
         col = list("#4DBBD5FF","#00A087FF"), 
         border = c("#2C738EFF","#4FC46AFF"), 
         main = "p-value = 3.985e-07",
         xlab = "", ylab = "Fraction of cells")
legend("topleft",bty="n",c("CD3+T cells", "Epithelial cells"),
       fill = c("#4DBBD5FF","#00A087FF"))
dev.off()