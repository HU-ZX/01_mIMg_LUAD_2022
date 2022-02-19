#### DEGs------
rm(list=ls())
library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
load("00.OriExp.Rda")
load("/Volumes/T7/Index_Gencode_v35/anno.Rdata")
rm(counts_df)
Targets<-data.table::fread("./Cluster.txt")
colnames(count0)<-substring(colnames(count0),1,12)

#### 1 Pre-processing----
### Organising sample information----
geneid<-data.frame(ENSG=rownames(count0))
genes<-merge(anno,geneid,by="ENSG")
genes<-genes[!duplicated(genes$Symbol),]
count0<-count0[genes$ENSG,]
### Organising gene annotations----
identical(colnames(count0),Targets$Sample)
x<-DGEList(count=count0,group=as.factor(Targets$Cluster))
group<-as.factor(Targets$Cluster)
x$samples$group
#x$samples$batch<-batch
x$samples
x$genes <- genes
x

#### 2 Data Processing ----
### 2.1 Transformations from the raw-scale----
L<-mean(x$samples$lib.size)* 1e-6
M<-median(x$samples$lib.size)* 1e-6
c(L,M)#56.21887 53.95298
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
summary(lcpm)


### 2.2 Removing genes that are lowly expressed----
table(rowSums(x$counts==0)==450)
# FALSE  TRUE 
# 26275   450
keep.exprs <- filterByExpr(x, group=as.factor(Targets$Cluster))
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)#  26275   450

pdf("./limma_edgeR_RNAseq123/gene_filter.pdf",height=10,width=10)
lcpm.cutoff <- log2(10/M + 2/L)
samplenames<-colnames(x$counts)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()

### 2.3 Normalising gene expression distributions----
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


### 2.4 Unsupervised clustering of samples----
lcpm <- cpm(x, log=TRUE)
pdf("./limma_edgeR_RNAseq123/pca.pdf",height=10,width=10)
par(mfrow=c(1,2))
col.group <- as.factor(Targets$Cluster)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, labels= as.factor(Targets$Cluster), col=col.group)

dev.off()

#### 3 Differential expression analysis----
### 3.1 Creating a design matrix and contrasts----
#design <- model.matrix(~0+group+batch)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
  C3vsC2= Cluster3 - Cluster2 , 
  C3vsC1= Cluster3 - Cluster1 ,
  C2vsC1= Cluster2 - Cluster1,
  levels = colnames(design))
contr.matrix
### 3.2 Removing heteroscedascity from count data----
pdf("./limma_edgeR_RNAseq123/Rm_heteroscedascity.pdf",height=5,width = 10)
par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
dev.off()
### 3.3 Examining the number of DE genes----
summary(decideTests(efit))
# C3vsC2 C3vsC1 C2vsC1
# Down     6118   4116   5756
# NotSig  12288  18985  15665
# Up       7869   3174   4854
dt <- decideTests(efit,lfc=0)
de.common <- which(dt[,1]!=0)
table(dt[,1]==1)
up<-which(dt[,1]==1)
table(dt[,1]==(-1))
down<-which(dt[,1]==(-1))
length(de.common)# 447
head(efit$genes$Symbol[de.common], n=20)
up<-efit[up,]
down<-efit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))
Df_tcga<-rbind(up,down)
Df_tcga <- Df_tcga[!is.na(Df_tcga$ENSG),] 
Df_tcga_pc<-Df_tcga%>%
  filter(Gene_type=="protein_coding")
all_tcga<- topTable(efit, coef=1, n=Inf,resort.by="logFC")
all_tcga<-all_tcga[!is.na(all_tcga$ENSG),]
save(Df_tcga,all_tcga,file = "./05Diffgene_tcga.rda")

up<-which(dt[,2]==1)
down<-which(dt[,2]==(-1))
up<-efit[up,]
down<-efit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))
Df_c3vsc1_tcga<-rbind(up,down)


up<-which(dt[,3]==1)
down<-which(dt[,3]==(-1))
up<-efit[up,]
down<-efit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))
Df_c2vsc1_tcga<-rbind(up,down)

Df_c3vsc2_tcga<-Df_tcga
Df_ls_tcga<-list(c3vsc2=Df_c3vsc2_tcga,
                 c3vsc1=Df_c3vsc1_tcga,
                 c2vsc1=Df_c2vsc1_tcga)
all_c3vsc1_tcga<- topTable(efit, coef=2, n=Inf,resort.by="logFC")
all_c2vsc1_tcga<- topTable(efit, coef=3, n=Inf,resort.by="logFC")
all_c3vsc2_tcga<- topTable(efit, coef=1, n=Inf,resort.by="logFC")
all_ls_tcga<-list(c3vsc2=all_c3vsc2_tcga,
                  c3vsc1=all_c3vsc1_tcga,
                  c2vsc1=all_c2vsc1_tcga)
save(Df_ls_tcga,all_ls_tcga,file="./05DiffGene_myDEG_tcga.rda")

### Volcano plot-----
library(ggrepel)
VolExp<-as.data.frame(all_tcga)
VolExp$color<-ifelse(VolExp$adj.P.Val<0.05 & abs(VolExp$logFC)>0,
                     ifelse(VolExp$logFC > 0,'red','blue'),
                     'black')

sub<-subset(VolExp,VolExp$Symbol%in%gene)
quantile(sub$logFC)
# 0%        25%        50%        75%       100% 
# -4.222929 -2.368769 -1.624377 -1.048837  0.286996 
table(sub$logFC>0.2)
table(VolExp$color)
# black  blue   red 
# 18428   158   126
color <- c(red = "red",black = "black",blue = "blue")
yMax=max(-log10(all_tcga$adj.P.Val))
yMax=ifelse(yMax>100,100,yMax)
xMax=max(abs(all_tcga$logFC))
xMax=ifelse(xMax>10,10,xMax)

ggplot(VolExp, aes(x=logFC, y=-log10(adj.P.Val),fill=color)) +
  geom_jitter(size=2,shape=21,color='black') +
  theme_classic()+
  scale_fill_manual(values = color)+
  labs(x="log2 (fold change)",y="-log10 (adj-p-value)") +
  #scale_y_continuous(limits = c(0,yMax+0.2),expand = c(0,0)) +
  #scale_x_continuous(limits = c(-1.5,xMax)) +
  #geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  #geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  geom_jitter(data=sub,aes(x=logFC, y=-log10(adj.P.Val)),size=3,shape=21,color='#d72323')+
  #geom_text(data=sub,aes(x=logFC+0.1,y=-log10(adj.P.Val)+0.1,label=symbol),size=3)+
  theme(legend.position = "none",
        axis.text=element_text(size=10),
        axis.line=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title=element_text(size=12,face="plain"),
        panel.background =element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank())
