library(GEOquery)
### Methylation-----
rm(list=ls())
pda_cluster<-data.table::fread("./Cluster.txt")
eset <- getGEO("GSE66836",destdir = './GSE66863/',AnnotGPL = T,getGPL = F)
beta.m <- exprs(eset[[1]])

pD.all <- pData(eset[[1]])
colnames(pD.all)
pD.all[1,]
files<-list.files("./GSE66863/GSE66836_RAW", full.names = T)
str_locate_all(files[1],"_")
substring(files[1],1,52)
Files<-data.frame(Files=files,
                  geo_accession=str_extract(files,"GSM\\d+"),
                  Array=sapply(strsplit(files,"_"),"[",4),
                  Slide=sapply(strsplit(files,"_"),"[",3),
                  Basename=substring(files,1,52))
Files<-Files%>%distinct(geo_accession,.keep_all = T)
pD.all <-pD.all[,c(1,2,8,10:15,23,24)]
colnames(pD.all)<-c("title","geo_accession",
                    "histology","Gender",
                    "Tissue","Stage_m",
                    "Tp53_status_m","EGFR_status_m","KRAS_status_m",
                    "Sentrix_ID","Sentrix_Position")
table(pD.all$histology)
pD.all<-pD.all%>%inner_join(Files,by="geo_accession")%>%
  filter(histology=="lung adenocarcinoma")%>%
  mutate(title=paste0("Tumor_",substring(str_extract(title,"Sample\\d+"),7,)))%>%
  inner_join(pda_cluster,by="title")
colnames(pD.all)
pD.all$title[1:10]
paste0("Tumor_",substring(str_extract(pD.all$title[1:5],"Sample\\d+"),7,))
pda_cluster$title[1]
pD<-pD.all[,c(1,2,12:24,26)]
colnames(pD)[c(2,3,7,10,15)]<-c("geo_accession_m","Files_m",
                                "geo_accession_e","Gender",
                                "Files_e")

beta.m=beta.m[,pD$geo_accession_m]
library(RColorBrewer)
library(limma)
library(minfi)
library(ChAMP)
densityPlot(beta.m, sampGroups=pD$Cluster, legend=FALSE)
legend("top", legend = levels(factor(pD$Cluster)), 
       text.col=brewer.pal(8,"Dark2"))
table(pD$Cluster)
library(minfi)
rgSet <- read.metharray.exp(targets=pD)
sampleNames(rgSet) <- pD$geo_accession_m
rgSet
detP <- detectionP(rgSet)
head(detP)
# examine mean detection p-values across all samples to identify any failed samples
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,1))
barplot(colMeans(detP), col=pal[factor(pD$Cluster)], las=2, 
        cex.names=0.8, ylim=c(0,0.06), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(pD$Cluster)), fill=pal,
       bg="white")
qcReport(rgSet, sampNames=Targets$Cluster, 
         pdf="./GSE66863_GSE66836_Data/GSE66836_qcReport.pdf")

### Filter-----
## minfi----
# remove poor quality samples
keep <- colMeans(detP) < 0.05
table(keep)

rgSet <- rgSet[,keep]
rgSet
# remove poor quality samples from targets data
pD <- pD[keep,]
# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)
mSetSq <- preprocessQuantile(rgSet) 
mSetRaw <- preprocessRaw(rgSet)

### Filtering
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
# FALSE   TRUE 
# 117848 367664 
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# romove sex  chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
# FALSE   TRUE 
# 8279 359385 
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# exclude cross reactive probes 
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt


densityPlot(rgSet, sampGroups=pD$Cluster,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(pD$Cluster)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSqFlt), sampGroups=pD$Cluster,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(pD$Cluster)),
       text.col=brewer.pal(8,"Dark2"))
save(rgSet,mSetSqFlt,mSetSq,pD,file="./GSE66863_GSE66836_Data/02Methylation_byminfi.Rda")



## ChAMP----
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)            
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)

library(ChAMP)
beta<-getBeta(rgSet)
mSetRaw <- preprocessRaw(rgSet)
mval<-getM(mSetRaw)
identical(colnames(beta),pD$geo_accession_m)
identical(colnames(detP),pD$geo_accession_m)

Myload<-champ.filter(beta=beta,
                     M=mval,
                     pd=pD,
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
sample_rm<-setdiff(pD$geo_accession_m,Myload$pd$geo_accession_m)
sample_rm<-pD[pD$geo_accession_m%in%sample_rm,]
probe<-rownames(Myload$beta)
# remove any probes that have failed in one or more samples
# keep <- rowSums(detP[rownames(Myload$beta),] < 0.01) == ncol(Myload$beta) 
# table(keep)
# keep <- !(probe %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
keep<-!(probe %in% xReactiveProbes$TargetID)
table(keep)
probe<-probe[keep]
beta_fl<-Myload$beta[probe,]
# na.names<-rownames(beta_fl[probe,])[apply(beta_fl[probe,],1,function(x) sum(is.na(x))>0)]
# beta_fl<-beta[probe,][!(rownames(beta[probe,])%in%na.names),]
beta_norm<-champ.norm(beta_fl,arraytype="450K",cores=5)
densityPlot(rgSet, sampGroups=pD$Cluster,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(pD$Cluster)),
       text.col=brewer.pal(8,"Dark2"))
identical(colnames(beta_norm),Myload$pd$geo_accession_m)
densityPlot(beta_norm, sampGroups=Myload$pd$Cluster,main="Beta values by CHAMP", legend=FALSE)
legend("top", legend = levels(factor(Myload$pd$Cluster)),
       text.col=brewer.pal(8,"Dark2"))
table(Myload$pd$Cluster)
pD<-Myload$pd
write.table(pD,"./Cluster.txt",sep="\t",col.names = T,row.names = F,quote=F)
save(beta_norm,Myload,pD,file="./GSE66863_GSE66836_Data/02Methylation_byChAMP.Rda")

pdf(file="./GSE66863_GSE66836_Data/Normalizeddensity.pdf",width=4,height=4)
densityPlot(beta_norm ,sampGroups = pD$Cluster,pal=c("#48466d","#3d84a8","#46cdcf"),
            main="GSE66836 Normalized Data", legend=FALSE,cex.main=1,)
legend("top", legend = levels(factor(pD$Cluster)), 
       text.col=c("#48466d","#3d84a8","#46cdcf"))
dev.off()
dat=t(beta_norm)
dat[1:4,1:4] 
library("FactoMineR")
library("factoextra")  
dat.pca <- PCA(dat , graph = FALSE) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = pD$Cluster, # color by groups
             palette = c("#48466d","#3d84a8","#46cdcf"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = NULL,
             title="GSE66836 PCA Plot"
)

ggsave('./GSE66836_PCA.pdf',height=4,width=4)
## DMP by ChAMP------
rm(list=ls())
load("./02Methylation_byChAMP.Rda")
library(ChAMP)
myQC<-champ.QC(Myload$beta,Myload$pd$Cluster)
myDMP <- champ.DMP(beta=beta_norm, pheno=pD$Cluster)
table(myDMP_c2vsc3$feature)
myDMP_c2vsc3 <- myDMP[[3]]

write.table(cbind(Probe=rownames(myDMP_c2vsc3),myDMP_c2vsc3),file="./GSE66836_myDMP.txt",sep="\t",quote=F,col.names = T,row.names = F)

myDMP_gse66836<-myDMP
# pCluster3-pCluster1
# pCluster1                  -1
# pCluster3                   1
# Contrasts
# Levels      pCluster2-pCluster1
# pCluster1                  -1
# pCluster2                   1
# Contrasts
# Levels      pCluster3-pCluster2
# pCluster2                  -1
# pCluster3                   1
save(myDMP_gse66836,file="./03DiffMeth_myDMP_gse66836.Rda")
load("./03DiffMeth_myDMP_gse66836.Rda")
myDMP<-myDMP_gse66836
DMP.GUI(DMP=myDMP[[3]],beta=beta_norm, pheno=pD$Cluster)

Diffgenes<-rownames(myDMP[[3]])[abs(myDMP[[3]]$logFC)>0]
set1<-intersect(Diffgenes,rownames(DMPs))
HDOM<-myDMP[[1]][myDMP[[1]]$deltaBeta>0,]  # Hydroxymethylation Analysis
myDMR <- champ.DMR(beta=beta_norm, pheno=pD$Cluster,method="Bumphunter",compare.group = c("Cluster2","Cluster3"))   
write.table(cbind(Probe=rownames(myDMR$BumphunterDMR),myDMR$BumphunterDMR),file="./GSE66836_MyDMR.txt",sep="\t",quote=F,col.names = T,row.names = F)

myGSEA <- champ.GSEA(beta=beta_norm,DMP=myDMP[[3]],
                     DMR=myDMR, arraytype="450K",adjPval=0.05, method="fisher")
myGSEA_66836<-myGSEA
save(myGSEA_66836,file="./03DiffMeth_myGSEA_66836.Rda")


## Figure 3 B Allmethy----
rm(list=ls())
load("./02Methylation_byChAMP.Rda")
methy<-apply(beta_norm,2,mean)
methy<-data.frame(geo_accession_m=names(methy),methy=methy) 
methy<-methy%>%
  inner_join(pD[,c("geo_accession_m","Cluster","Tp53_status","EGFR_status","KRAS_status")],by="geo_accession_m") %>% 
  arrange("Cluster")
my_comparisons <- list(c("Cluster1", "Cluster2"), c("Cluster1", "Cluster3"), c("Cluster2", "Cluster3"))
with(methy,kruskal.test(methy~Cluster))
with(methy,tapply(methy,Cluster,shapiro.test))
bartlett.test(methy~Cluster,data = methy) 
methy_kras<-methy %>% filter(!(is.na(Tp53_status)))
pdf("./Cluster_Phenotype_Plot/02Methylation_bytp53.pdf",width=6,height=4)
ggplot(methy_kras,aes(Cluster,methy,fill=Cluster))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap("Tp53_status",ncol=2)+
  labs(x = "", y = "Average methylation level")+
  theme(axis.text.x = element_text(colour = "black",size=12),
        axis.text.y = element_text(colour = "black",size=10),
        axis.title.y = element_text(colour = "black",size=12),
        strip.background=element_blank(),
        strip.text = element_text(colour = "black",size=12))  +
  theme(axis.ticks.x = element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c('#48466d','#3d84a8','#46cdcf'),labels=c('Cluster 1','Cluster 2','Cluster 3'))+
  scale_y_continuous()+ 
  stat_compare_means(comparisons=my_comparisons,label="p.signif")+  
  stat_compare_means(label.y=0.65)

dev.off()
## Figure S1 B PCA plot----
load("./01Importexpression.Rda")
load("./02Methylation_byChAMP.Rda")
Target<-data.table::fread("./Cluster.txt")
pD<-Target %>% 
  column_to_rownames("geo_accession_m")
pD<-pD[colnames(beta_norm),]

exp<-expr_norm[,pD$geo_accession_e]


load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"gene"])
dat1<-t(exp[probeset1$gene,])
dat2<-t(beta_norm[probeset1$Probe,])
dat[1:4,1:4] 
library("FactoMineR")
library("factoextra")  
dat.pca <- PCA(dat2 , graph = FALSE) 
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

