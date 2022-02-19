####  Construct RGSet-----
Targets<-data.table::fread("./Cluster.txt")
library(minfi)
rgSet <- read.metharray.exp(targets=Targets)
sampleNames(rgSet) <- Targets$Sample
rgSet
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
# examine mean detection p-values across all samples to identify any failed samples
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,1))
barplot(colMeans(detP), col=pal[factor(Targets$Cluster)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(Targets$Cluster)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(Targets$Cluster)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(Targets$Cluster)), fill=pal, 
       bg="white")
qcReport(rgSet, sampNames=Targets$Cluster, 
         pdf="qcReport.pdf")
# remove poor quality samples
keep <- colMeans(detP) < 0.05
table(keep)
rgSet <- rgSet[,keep]
rgSet
# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)
save(rgSet,detP,Targets,file="01.TCGA_Meth_Clinic.Rda")



### CHAMP line----
rm(list=ls())
load("01.TCGA_Meth_Clinic.Rda")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k_df<-ann450k%>%as.data.frame()
head(ann450k)
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
library(ChAMP)
beta<-getBeta(rgSet)
mSetRaw <- preprocessRaw(rgSet)
mval<-getM(mSetRaw)
identical(colnames(beta),Targets$Sample)
identical(colnames(detP),Targets$Sample)
Myload<-champ.filter(beta=beta,
                     M=mval,
                     pd=Targets,
                     detP=detP,
                     autoimpute=TRUE,
                     filterDetP=TRUE,
                     ProbeCutoff=0,
                     SampleCutoff=0.1,
                     detPcut=0.05,
                     filterBeads=TRUE,
                     filterNoCG = TRUE,
                     filterSNPs = TRUE,
                     filterMultiHit = TRUE,
                     filterXY = TRUE,
                     fixOutlier = TRUE,
                     arraytype = "450K")
str(Myload)
probe<-rownames(Myload$beta)
# remove any probes that have failed in one or more samples
# keep <- rowSums(detP[rownames(Myload$beta),] < 0.01) == ncol(Myload$beta) 
# table(keep)
# keep <- !(probe %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
keep<-!(probe %in% xReactiveProbes$TargetID)
table(keep)
probe<-probe[keep]
beta_fl <- Myload$beta[probe,]
beta_norm<-champ.norm(beta_fl,arraytype="450K",cores=5)
save(beta_norm,Targets,Myload,file="04CHAMPline_import.Rda")

pdf(file="./Normalizeddensity.pdf",width=4,height=4)
densityPlot(beta_norm ,sampGroups=Targets$Cluster,pal=c("#48466d","#3d84a8","#46cdcf"),
            main="TCGA LUAD Normalized Data", legend=FALSE,cex.main=1,)
legend("top", legend = levels(factor(Targets$Cluster)), 
       text.col=c("#48466d","#3d84a8","#46cdcf"))
dev.off()

# densityPlot(beta_norm ,sampGroups=Targets$Cluster, 
#             main="Normalized Beta values by ChAMP", legend=FALSE)
# legend("top", legend = levels(factor(Targets$Cluster)), 
#        text.col=brewer.pal(8,"Dark2"))
# densityPlot(bVals, sampGroups=Targets$Cluster, main="Normalized Beta values by minfi", 
#             legend=FALSE, xlab="Beta values")
# legend("top", legend = levels(factor(Targets$Cluster)), 
#        text.col=brewer.pal(8,"Dark2"))

dat=t(beta_norm)
dat[1:4,1:4] 
library("FactoMineR")
library("factoextra")  
dat.pca <- PCA(dat , graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Targets$Cluster, # color by groups
             palette = c("#48466d","#3d84a8","#46cdcf"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = NULL,
             title="TCGA LUAD PCA Plot"
)
ggsave('./TCGALUAD_PCA.pdf',height=4,width=4)
###DMP by ChAMP----
rm(list=ls())
load("04CHAMPline_import.Rda")
library(ChAMP)
myDMP <- champ.DMP(beta=beta_norm, pheno=Targets$Cluster)
table(myDMP_c2vsc3$feature)
myDMP_c2vsc3 <- champ.DMP(beta=beta_norm, pheno=Targets$Cluster,compare.group=c("Cluster2","Cluster3"))[[1]]
DMP.GUI(DMP=myDMP[[3]],beta=beta_norm, pheno=Targets$Cluster)
write.table(cbind(Probe=rownames(myDMP_c2vsc3),myDMP_c2vsc3),file="./TCGA_MyDMP.txt",sep="\t",quote=F,col.names = T,row.names = F)
myDMP_tcga<-myDMP
# Contrasts
# Levels      pCluster2-pCluster1
# pCluster1                  -1
# pCluster2                   1
# Contrasts
# Levels      pCluster3-pCluster1
# pCluster1                  -1
# pCluster3                   1
# Contrasts
# Levels      pCluster3-pCluster2
# pCluster2                  -1
# pCluster3                   1
save(myDMP_tcga,file="./03DiffMeth_myDMP_tcga.Rda")



Diffgenes<-rownames(myDMP[[3]])[abs(myDMP[[3]]$logFC)>0.2]
set1<-intersect(Diffgenes,rownames(DMPs))
HDOM<-myDMP[[1]][myDMP[[1]]$deltaBeta>0,]  # Hydroxymethylation Analysis
myDMR <- champ.DMR(beta=beta_norm, pheno=Targets$Cluster,compare.group=c("Cluster2","Cluster3"),method="Bumphunter")   
write.table(cbind(Probe=rownames(myDMR$BumphunterDMR),myDMR$BumphunterDMR),file="./TCGA_MyDMR.txt",sep="\t",quote=F,col.names = T,row.names = F)

myGSEA <- champ.GSEA(beta=beta_norm,DMP=myDMP[[3]],
                     DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
myGSEA_tcga<-myGSEA
save(myGSEA_tcga,file="./03DiffMeth_myGSEA_tcga.Rda")
require(Biobase)
deg<-myDMP[[3]]
head(deg)
deg$g<-ifelse(abs(deg$logFC) < 0.2,'stable',
              ifelse(deg$logFC > 0.2,'UP','DOWN'))
table(deg$g)
deg$symbol<-deg$gene






### Figure 3 A All methylation----
methy<-apply(beta_norm,2,mean)
methy<-data.frame(Sample=names(methy),methy=methy) 
methy<-methy%>%
  inner_join(target[,c("Sample","Cluster","egfr_status")],by="Sample") %>% 
  arrange("Cluster")
with(methy,tapply(methy,Cluster,shapiro.test))
bartlett.test(methy~Cluster,data = methy) 
wilcox.test(methy~KRAS_status,methy)
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))

pdf("./Cluster_Phenotype_Plot/02Methylation.pdf",width=6,height=4)
ggplot(methy,aes(Cluster,methy,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  #facet_wrap("egfr_status",ncol=2)+
  labs(x = "", y = "Average methylation level")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12),
        strip.background=element_blank(),
        strip.text = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'),labels=c('Cluster 1','Cluster 2','Cluster 3'))+
  scale_y_continuous()+ 
  stat_compare_means(comparisons=my_comparisons,label="p.signif")+  
  stat_compare_means(label.y=0.65)
dev.off()

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



library(GSVA)
library(tidyverse)
library(ggpubr)
Sig<-data.table::fread("/Volumes/T7/Methylation/02Signature_gene&probe.txt") %>% 
  arrange(gene)
Sig<-data.table::fread("/Volumes/T7/Methylation/02tcga_Signature_gene&probe.txt") %>% 
  arrange(gene)
gene<-unique(Sig$gene)
### Figure S1 PCAplot-----

load("02.TPM&Targets.Rda")
load("04CHAMPline_import.Rda")
Target<-data.table::fread("./Cluster_mut.txt")
load("/Users/huzixin/Data_Project/Methylation/00Gene_exp_methy_match.rda")
probe_exp_meth<-rbind(probe_exp_methy_d_cis,
                      probe_exp_methy_d_trans,
                      probe_exp_methy_u_cis,
                      probe_exp_methy_u_trans
                      
)
load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
dat=t(beta_norm[probe_exp_meth$Probe,])
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"gene"])

dat1<-t(beta_norm[probeset1$Probe,])
dat2<-t(TPM_sym[probeset1$gene,])
dat[1:4,1:4] 
library("FactoMineR")
library("factoextra")  
dat.pca1 <- PCA(dat1 , graph = FALSE) 
dat.pca2 <- PCA(dat2 , graph = FALSE) 
pdf("./Cluster_Phenotype_Plot/methy_exp_pca.pdf",width=4,height=4)
fviz_pca_ind(dat.pca1,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Targets$Cluster, # color by groups
             palette = c("#48466d","#3d84a8","#46cdcf"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = NULL)+
  theme(axis.text = element_text(colour = "black",size=12),
        axis.title = element_text(colour = "black",size=12)) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.title = element_blank())
dev.off()
