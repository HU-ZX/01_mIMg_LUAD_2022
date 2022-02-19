# GSE60644 GPL10558
### Expression data----
library(tidyverse)
require("illuminaHumanv4.db")
library("beadarray")
library("GEOquery")
library("limma")

files<-list.files("./GSE60645/GSE60644/", full.names = T)
target<-getGEO(GEO="GSE60644",file=files[8],getGPL=F)
pda<-pData(target)
files_sup<-list.files("./GSE60645/GSE94601/", full.names = T)
target_sup<-getGEO(GEO="GSE94601",file=files_sup[3],getGPL=F)
pda_sup<-pData(target_sup)
colnames(pda_sup)
colnames(pda)
pda_co<-pda%>%
  dplyr::select("title","geo_accession")%>%
  mutate(title=str_extract(pda$title,"LU\\d+"))%>%
  inner_join(pda_sup[,c(1,2,10:12)],by="title")
colnames(pda_co)[4:6]<-c("age","histology","sex")
pda_co<-pda_co%>%
  mutate(age=as.numeric(sapply(strsplit(as.character(age),":"),"[",2)),
         histology=substring(sapply(strsplit(as.character(histology),":"),"[",2),2,),
         sex=substring(sapply(strsplit(as.character(sex),":"),"[",2),2,))  
non_norm<- read.ilmn(files = files[6],
                                 probeid = "ID_REF")
head(non_norm$E)
dim(non_norm$E)
boxplot(log2(non_norm$E),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

library(arrayQualityMetrics)

rownames(pda_co) <-  pda_co$title
identical(colnames(non_norm$E),pda_co$title)
eset <- ExpressionSet(assayData = non_norm$E[,pda_co$title], 
                      phenoData = AnnotatedDataFrame(data = pda_co))

arrayQualityMetrics(eset, 
                    outdir = "./GSE60645/GSE60644/GSE60644_QC_aqm", 
                    do.logtransform = TRUE,
                    force = TRUE,
                    intgroup = "histology")
dev.off()
norm<-neqc(non_norm, detection.p="Detection")
norm<-norm$E[,pda_co$title]
boxplot(log2(norm),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)

x <- illuminaHumanv4SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
xxx<-sapply(xx,as.character)%>%as.data.frame()%>%
  rownames_to_column("ID_REF")
colnames(xxx)[2]<-"Gene_symbol"
table(pda_co$histology)
head(xxx)
PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  object <- expr
  rv <- genefilter::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
                  PC2 = pca$x[, 2], 
                  group = group, 
                  name = colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
                      size = 3,
                      segment.color = "black",
                      show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}
PCA_new(norm, 
        ntop = nrow(norm),
        group = pda_co$histology)
dim(norm)
norm<-norm%>%as.data.frame()%>%
  rownames_to_column("ID_REF")%>%
  inner_join(xxx,by="ID_REF")%>%dplyr::select(-"ID_REF")%>%
  aggregate(.~Gene_symbol,data=.,max)%>%
  column_to_rownames("Gene_symbol")

save(file="./GSE60644_GSE60645_Data/01Importdata.rda",pda_co,non_norm,norm,xxx)
