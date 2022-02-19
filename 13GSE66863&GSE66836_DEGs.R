### DEGs-----
rm(list=ls())
library(tidyverse)

pda<-data.table::fread("./Cluster.txt")
load("./01Importexpression.Rda")
expr_norm<-expr_norm %>% column_to_rownames("Gene")
expr_norm<-expr_norm[,pda$geo_accession]

library(limma)
group <- pda$Cluster
design <- model.matrix(~ 0 + group, data = pda)
colnames(design) <- gsub("group", "", colnames(design))
head(design)

## contrasts
contr.matrix <- makeContrasts(
  C3vsC2= Cluster3 - Cluster2 , 
  C3vsC1= Cluster3 - Cluster1 ,
  C2vsC1= Cluster2 - Cluster1,
  levels = colnames(design))
contr.matrix

fit <- lmFit(expr_norm, design)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit)
summary(decideTests(fit))
# C3vsC2 C3vsC1 C2vsC1
# Down     3494   2338    599
# NotSig  23209  27280  32631
# Up       8049   5134   1522
dt <- decideTests(fit,lfc=0)
de.common <- which(dt[,1]!=0)
table(dt[,1]==1)
up<-which(dt[,1]==1)
table(dt[,1]==(-1))
down<-which(dt[,1]==(-1))
length(de.common)# 447

up<-fit[up,]
down<-fit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))

Df_gse66863<-rbind(up,down)%>%as.data.frame()%>%rownames_to_column("Symbol")
table(Df_gse66863$up_down)
all_gse66863<-topTable(fit, coef=1, n=Inf,resort.by="logFC")
save(Df_gse66863,all_gse66863,file = "./05Diffgene_gse66863.rda")


up<-which(dt[,2]==1)
down<-which(dt[,2]==(-1))
up<-fit[up,]
down<-fit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))
Df_c3vsc1_gse66863<-rbind(up,down)
table(Df_c3vsc1_gse66863$up_down)

up<-which(dt[,3]==1)
down<-which(dt[,3]==(-1))
up<-fit[up,]
down<-fit[down,]
up <- topTable(up, coef=1, n=Inf,resort.by="logFC")
up$up_down<-rep("up",nrow(up))
down <- topTable(down, coef=1, n=Inf,resort.by="logFC")
down$up_down<-rep("down",nrow(down))
Df_c2vsc1_gse66863<-rbind(up,down)
table(Df_c2vsc1_gse66863$up_down)

Df_c3vsc2_gse66863<-Df_gse66863
Df_ls_gse66863<-list(c3vsc2=Df_c3vsc2_gse66863,
                     c3vsc1=Df_c3vsc1_gse66863,
                     c2vsc1=Df_c2vsc1_gse66863)
all_c3vsc1_gse66863<- topTable(fit, coef=2, n=Inf,resort.by="logFC")
all_c2vsc1_gse66863<- topTable(fit, coef=3, n=Inf,resort.by="logFC")
all_c3vsc2_gse66863<- topTable(fit, coef=1, n=Inf,resort.by="logFC")
all_ls_gse66863<-list(c3vsc2=all_c3vsc2_gse66863,
                      c3vsc1=all_c3vsc1_gse66863,
                      c2vsc1=all_c2vsc1_gse66863)
save(Df_ls_gse66863,all_ls_gse66863,file="./05DiffGene_myDEG_gse66863.rda")
