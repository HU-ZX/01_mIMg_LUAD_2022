##Import Methylation----
rm(list=ls())
library(tidyverse)
library("GEOquery")
library("limma")
library("Biobase")
library(ChAMP)
library("minfi")
library(RColorBrewer)
pda_cluster<-data.table::fread("./Cluster.txt")
eset <- getGEO("GSE56044",destdir = './GSE60645/GSE56044/', AnnotGPL = T,getGPL = F)

beta.m <- exprs(eset[[1]])
pD.all <- pData(eset[[1]])
pD.all[1,]
colnames(pD.all)
pD <- pD.all[, c("title", "geo_accession")]
colnames(pD)[1]<-"Sample"
pD <- pD%>%
  inner_join(pda_cluster,by="Sample")%>%
  arrange(Cluster)
rownames(pD)<-pD$title

beta.m <-beta.m [,pD$geo_accession]

colnames(beta.m)<-pD$title
dim(beta.m)
na.names<-rownames(beta.m)[apply(beta.m,1,function(x) sum(is.na(x))>0)]
beta.m<-beta.m[!(rownames(beta.m)%in%na.names),]
table(is.na(beta.m))
densityPlot(beta.m,sampGroups = pD$Cluster,
            main="Normalized", legend=FALSE)

head(pD)
myLoad<- champ.filter(beta = beta.m ,pd = pD)
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
probe<-rownames(myLoad$beta)
# remove any probes that have failed in one or more samples
# keep <- rowSums(detP[rownames(Myload$beta),] < 0.01) == ncol(Myload$beta) 
# table(keep)
# keep <- !(probe %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
keep<-!(probe %in% xReactiveProbes$TargetID)
table(keep)
probe<-probe[keep]
beta_fl <- myLoad$beta[probe,]
myLoad$beta<-beta_fl
CpG.GUI(CpG=rownames(myLoad$beta),arraytype="450K")
save(myLoad,file="./02myLoad.Rda")

## CHAMP Normlization Line----
load("./02myLoad.Rda")
myQC<-champ.QC(myLoad$beta,myLoad$pd$Cluster)
densityPlot(myLoad$beta ,sampGroups = pD$Cluster,pal=c("#48466d","#3d84a8","#46cdcf"),
            main="GSE56044 Normalized Data", legend=FALSE,cex.main=1)
QC.GUI(beta=myLoad$beta,myLoad$pd$Cluster,arraytype="450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
mySVD<-champ.SVD()
?densityPlot()
#c(Cluster1 ="#48466d",Cluster2="#3d84a8",Cluster3="#46cdcf")
pdf(file="./GSE60644_GSE60645_Data/Normalizeddensity.pdf",width=4,height=4)
densityPlot(myNorm ,sampGroups = pD$Cluster,pal=c("#48466d","#3d84a8","#46cdcf"),
            main="GSE56044 Normalized Data", legend=FALSE,cex.main=1)
legend("top", legend = levels(factor(pD$Cluster)), 
       text.col=c("#48466d","#3d84a8","#46cdcf"))
dev.off()

save(myNorm,pD,file= "./02Methylation.Rda")

### Figure S1 C PCA----
load("01Importdata.rda")
load("02Methylation.Rda")
Target<-data.table::fread("./Cluster.txt")
beta<-myNorm[,Target$Sample]
exp<-norm[,Target$Sample]
load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"gene"])
dat=t(beta[probeset1$Probe,])
dat[1:4,1:4] 
library("FactoMineR")
library("factoextra")  
dat.pca <- PCA(dat , graph = FALSE) 
pdf("./Cluster_Phenotype_Plot/methy_exp_pca.pdf",width=4,height=4)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Target$Cluster, # color by groups
             palette = c("#48466d","#3d84a8","#46cdcf"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = NULL)+
  theme(axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()


## DMPs-----
library(ChAMP)
load("02Methylation.Rda")
table(pD$Cluster)
myDMP <- champ.DMP(beta=myNorm, pheno=pD$Cluster)
myDMP_c2vsc3<-myDMP[[2]]
table(myDMP_c2vsc3$logFC<0)
myDMP_c2vsc3<-champ.DMP(beta=myNorm, pheno=pD$Cluster,compare.group=c("Cluster2","Cluster3"))
write.table(cbind(Probe=rownames(myDMP_c2vsc3),myDMP_c2vsc3),file="./GSE56044_MyDMP.txt",sep="\t",quote=F,col.names = T,row.names = F)
myDMP_gse56044<-myDMP
save(myDMP_gse56044,file="./03DiffMeth_myDMP_gse56044.Rda")

head(myDMP[[1]])
dim(myDMP[[1]])
myDMP[[1]][myDMP[[1]]$deltaBeta>0,]  # Hydroxymethylation Analysis
DMP.GUI(DMP=myDMP[[1]],beta=beta.m, pheno=pD$im_f)

### Figure 3 C All methy----
load("02Methylation.Rda")
methy<-apply(myNorm,2,mean)
methy<-data.frame(title=names(methy),methy=methy) 
methy<-methy%>%
  inner_join(pD[,c("title","Cluster")],by="title") %>% 
  arrange("Cluster")
library(tidyverse)
library(ggpubr)
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
with(methy,tapply(methy,Cluster,shapiro.test))
bartlett.test(methy~Cluster,data = methy) 
pdf("./Cluster_Phenotype_Plot/02Methylation_all.pdf",width=4,height=4)
ggplot(methy,aes(Cluster,methy,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  #facet_grid("Cluster")+
  labs(x = "", y = "Average methylation level")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'),labels=c('Cluster 1','Cluster 2','Cluster 3'))+
  scale_y_continuous()+ 
  stat_compare_means(comparisons=my_comparisons,label="p.signif")+  
  stat_compare_means()

dev.off()



