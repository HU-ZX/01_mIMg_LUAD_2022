### Figure 4 A Enrichment------
rm(list=ls())
load("./00Gene_exp_methy_match.rda")
gene1<-unique(probe_exp_methy_d_cis$gene)
gene2<-unique(probe_exp_methy_d_trans$gene)
gene3<-unique(probe_exp_methy_u_cis$gene)
gene4<-unique(probe_exp_methy_u_trans$gene)

gene_1<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_u_trans$gene))
gene_2<-unique(c(probe_exp_methy_d_trans$gene,probe_exp_methy_u_cis$gene))
gene_all<-unique(probe_exp_meth$gene)
gene_1<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_d_trans$gene))
gene_2<-unique(c(probe_exp_methy_u_trans$gene,probe_exp_methy_u_cis$gene))
library(clusterProfiler)
library(org.Hs.eg.db)
gene.df <- bitr(gene_all,
                fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)# 
ego_CC<- enrichGO(gene = gene.df$ENTREZID,
                  OrgDb = org.Hs.eg.db, 
                  ont   = "CC",
                  pvalueCutoff =0.05, 
                  qvalueCutoff =0.05,
                  readable = TRUE)
ego_BP@result$Cat<-"BP"
save(ego_1_BP,ego_1_CC,ego_1_MF,ego_2_BP,ego_2_CC,ego_2_MF,file="./DiffDMP&DEG_enrichment2.rda")
save(ego_BP,ego_CC,ego_MF,file="./DiffDMP&DEG_enrichment3.rda")



a<-filter(arrange(ego_BP@result,desc(Count),pvalue),pvalue<0.01)
b<-filter(arrange(ego_MF@result,desc(Count),pvalue),pvalue<0.01)
c<-filter(arrange(ego_CC@result,desc(Count),pvalue),pvalue<0.01)
write.csv(C,"./TableS3_1CC.csv")


rm(list=ls())
load("./DiffDMP&DEG_enrichment3.rda")
# result<-rbind(arrange(ego_1_BP@result,desc(Count),pvalue)[1:6,],
#               arrange(ego_1_MF@result,desc(Count),pvalue)[1:6,],
#               arrange(ego_1_CC@result,desc(Count),pvalue)[1:6,]) 
# 
# result<-rbind(
#               arrange(ego_2_BP@result,desc(Count),pvalue)[1:6,],
#               arrange(ego_2_MF@result,desc(Count),pvalue)[1:6,],
#               arrange(ego_2_CC@result,desc(Count),pvalue)[1:6,])

result<-rbind(
  arrange(ego_BP@result,desc(Count),pvalue)[1:10,],
  arrange(ego_MF@result,desc(Count),pvalue)[1:10,],
  arrange(ego_CC@result,desc(Count),pvalue)[1:10,])  



gene<-strsplit(result$geneID,"/")

count<-sapply(1:length(gene),function(x) length(gene[[x]]))
genelist<-probe_exp_meth[,c("gene","logFC_g_tcga")] %>% 
  distinct(gene,.keep_all = T)
logFC<-sapply(unlist(gene),function(x) genelist$logFC_g_tcga[match(x,genelist$gene)])
s<-1
zsc<-c()
for(c in 1:length(count)){
  value<-0
  e<-s+count[c]-1
  value<-sapply(logFC[s:e],function(x) ifelse(x>0,1,-1))
  zsc<-c(zsc,sum(value)/sqrt(count[c]))
  s<-e+1
}
result$Scaled_enrichment_score<-zsc
result$Status<-ifelse(result$Scaled_enrichment_score>0,"UP","DOWN")

result<-result%>%
  arrange(Cat,desc(abs(Scaled_enrichment_score)),desc(Count),desc(pvalue))
a<-setdiff(result1$ID,result$ID)
a<-result1[result1$ID%in%a,]
a<-a %>% arrange(desc(abs(Scaled_enrichment_score)),desc(Count),desc(pvalue))
library(scales)
library(RColorBrewer)
mypal<-pal_lancet("lanonc")(8)
brewer.pal(5,"Pastel1")
mypal=c('#df5e88','#f6efa6','#fce2ce','#a2d5f2','#7ac7c4')
show_col(mypal)
print(mypal)

pdf("./Phenotype/all.pdf",height=8,width=8)
ggplot(result, aes(Scaled_enrichment_score, Description,fill=Count,size=-log10(pvalue))) + 
  facet_wrap('Cat',scales = "free",nrow=3)+
  theme_minimal()+
  labs(x = "Enrichment score", y = "") + 
  geom_point(shape = 21, col = "black", alpha = 0.5) + 
  scale_fill_continuous(type = "viridis")+
  theme(axis.text=element_text(size=8),
        axis.line=element_line(color="grey80"),
        axis.ticks=element_line(color="grey80"),
        axis.title=element_text(size=10,face="bold"),
        panel.background =element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        strip.text = element_text(size=8)
  )+
  guides(fill=guide_colorbar("Count",order=1))
dev.off()


## Figure 4 C Circos plot of corrrelation ----
library(OmicCircos)
library(tidyverse)
library(scales)
show_col(c('#F8766D','#7cae00','#01bfc4','#c77cff'))
#'#F8766D' cis hyper
#'#7cae00' cis hypo
#'#01bfc4' trans hyper
#'#c77cff' trans hypo

rm(list=ls())
load("./00Probe_filtercr_sp.rda")
probe_exp_meth_cr$cr[probe_exp_meth_cr$logFC_m_tcga>0&probe_exp_meth_cr$logFC_g_tcga>0]<-'#F8766D'
probe_exp_meth_cr$cr[probe_exp_meth_cr$logFC_m_tcga<0&probe_exp_meth_cr$logFC_g_tcga<0]<-'#7cae00'
probe_exp_meth_cr$cr[probe_exp_meth_cr$logFC_m_tcga>0&probe_exp_meth_cr$logFC_g_tcga<0]<-'#01bfc4'
probe_exp_meth_cr$cr[probe_exp_meth_cr$logFC_m_tcga<0&probe_exp_meth_cr$logFC_g_tcga>0]<-'#c77cff'

Gene.Label<-probe_exp_meth_cr[,c("Probe","gene","Start","End","Chromosome",
                                 "cr","Estimate_tcga","Estimate_gse56044","Estimate_gse66836")]%>%
  mutate(Chromosome=substring(as.character(Chromosome),4),
         po=Start) 
chr_order<-as.character(unique(str_sort(as.numeric(Gene.Label$Chromosome),numeric = T)))
b<-data.frame()
for(i in 1:length(chr_order)){
  a<-Gene.Label[Gene.Label$Chromosome==chr_order[i],]%>%arrange(po)
  b<-rbind(b,a)
}
Gene.Label<-b
colnames(Gene.Label)[5]<-"chr"
rm(a,b)

map1<-Gene.Label[,c(5,10,2)]%>%distinct(gene,.keep_all = T)
colnames(map1)<-c("chr","po","Gene")
map2<-Gene.Label[,c(5,10,1)]%>%distinct(Probe,.keep_all = T)
colnames(map2)<-c("chr","po","Probe");
Gene.Label$es.color<-ifelse(Gene.Label$Estimate_tcga>0,"red","blue")

es1<-Gene.Label[,c("chr","po","Probe","Estimate_tcga","es.color")]
es2<-Gene.Label[,c("chr","po","Probe","Estimate_gse56044","es.color")]
es3<-Gene.Label[,c("chr","po","Probe","Estimate_gse66836","es.color")]

#tcga----
load("./Download/02.TPM&Targets.Rda")
Targets<-data.table::fread("./Download/Cluster.txt")%>%arrange(Cluster)
exp1<-log2(TPM_sym[map1$Gene,Targets$Sample]+1)%>%
  as.data.frame()
exp1<-cbind(map1,exp1)
colnames(exp1)[3]<-"NAME"
load("./Download/04CHAMPline_import.Rda")
meth1<-beta_norm[Gene.Label$Probe,Targets$Sample]%>%
  as.data.frame()
meth1<-cbind(map2,meth1)
colnames(meth1)[3]<-"NAME"
rm(beta_norm,Myload,Targets)
rm(TPM,TPM_sym,Targets)
#gse66863----
load("./GEO/GSE66863_GSE66836_Data/01Importexpression.Rda")
load("./GEO/GSE66863_GSE66836_Data/02Methylation_byChAMP.Rda")
Targets<-data.table::fread("./GEO/GSE66863_GSE66836_Data/Cluster.txt")%>%arrange(Cluster)
exp2<-expr_norm%>%column_to_rownames("Gene")
exp2<-exp2[map1$Gene,Targets$geo_accession_e]
exp2<-cbind(map1,exp2)
colnames(exp2)[3]<-"NAME"
rm(expr_norm,norm,pda,GSE66863_raw)

meth2<-beta_norm[map2$Probe,Targets$geo_accession_m]%>%
  as.data.frame()
meth2<-cbind(map2,meth2)
colnames(meth2)[3]<-"NAME"
rm(beta_norm,Myload,pD)
rm(Targets)

#gse60644----
load("./GEO/GSE60644_GSE60645_Data/02Methylation.Rda")
load("./GEO/GSE60644_GSE60645_Data/01Importdata.rda")
Targets<-pD
exp3<-norm[map1$Gene,Targets$title]
exp3<-cbind(map1,exp3)
colnames(exp3)[3]<-"NAME"
rm(non_norm,norm,pda_co,xxx)

meth3<-myNorm[map2$Probe,Targets$title]%>%
  as.data.frame()
meth3<-cbind(map2,meth3)
colnames(meth3)[3]<-"NAME"
rm(myNorm,pD)
rm(Targets)

# Plot----
library(OmicCircos)
pdf("./Methy_Diff/01probe_exp_meth_cr_Circos.pdf",width=8,height=8)
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
circos(R=360, type="chr", cir="hg19", print.chr.lab=T, W=4);
circos(R=370, cir="hg19", W=20, mapping=map1, type="label", 
       side="out", col="black",cex=0.6);
colors<-Gene.Label$cr;
circos(R=350, cir="hg19", W=20, mapping=map2, type="label", 
       side="in",cex=0.6,col=colors);

# bk<-c(seq(-0.9,0,by = 0.01),seq(0+0.01,0.9,by=0.01))
# col_meth=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
#            colorRampPalette(colors=c("white","red"))(length(bk)/2))

circos(R=240, cir="hg19", W=30, 
       mapping=es1, cex=1,
       col.v=4, type="s",col=es1$es.color, B=T, scale=TRUE)


#circos(R=200, cir="hg19", W=20, mapping=exp1,type="heatmap2",cluster = F,col.bar =T, lwd=0.1, col="blue");
#circos(R=220, cir="hg19", W=30, mapping=meth1,type="heatmap2",cluster = F,col.bar =F, lwd=0.1);


circos(R=200, cir="hg19", W=30, 
       mapping=es2, cex=1,
       col.v=4, type="s",col=es2$es.color, B=T, scale=TRUE)

#circos(R=200, cir="hg19", W=30, mapping=exp3,type="heatmap2",cluster = F,col.bar =F, lwd=0.1, col="blue");
#circos(R=140, cir="hg19", W=30, mapping=meth3,type="heatmap2",cluster = F,col.bar =F, lwd=0.1, col="blue");

circos(R=160, cir="hg19", W=30, 
       mapping=es3, cex=1,
       col.v=4, type="s",col=es3$es.color, B=T, scale=TRUE)
#circos(R=140, cir="hg19", W=20, mapping=exp2,type="heatmap2",cluster = F,col.bar =T, lwd=0.1, col="blue");
#circos(R=60, cir="hg19", W=30, mapping=meth2,type="heatmap2",cluster = F,col.bar =F, lwd=0.1, col="blue");
dev.off()



## plot Circos  ----
rm(list=ls())
load("/Users/huzixin/Data_Project/Methylation/01IMps_Pheno_Probe_sp.rda")
load("/Users/huzixin/Data_Project/Methylation/00Probe_filtercr_sp.rda")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%probeset1$Probe,"gene"])
colnames(probe_exp_meth_cr)
probe_exp_meth_cr<-probe_exp_meth_cr[,c("Probe","gene","Chromosome","Start","End",
                                        "feature","cgi","Enhancer","UCSC_CpG_Islands_Name",
                                        "logFC_m_tcga","logFC_g_tcga","Estimate_tcga",
                                        "logFC_m_gse60644","logFC_g_gse60644","Estimate_gse56044",
                                        "logFC_m_gse66836","logFC_g_gse66836","Estimate_gse66836")] %>% 
  arrange(desc(gene))
table(probe_exp_meth_cr$cgi)

write.table(probe_exp_meth_cr,file="./TableS3probe_exp_meth_cr.txt",sep="\t",quote=F,col.names = T,row.names=F)
write.table(probeset1,file="./TableSmIMg24_21.txt",sep="\t",quote=F,col.names = T,row.names=F)
Gene.Label<-probe_exp_meth_cr[,c("Probe","gene","Start","End","Chromosome",
                                 "cr","Estimate_tcga","Estimate_gse56044","Estimate_gse66836")]%>%
  mutate(Chromosome=substring(as.character(Chromosome),4),
         po=Start) %>% 
  filter(Probe%in%probeset1$Probe)
chr_order<-as.character(unique(str_sort(as.numeric(Gene.Label$Chromosome),numeric = T)))
b<-data.frame()
for(i in 1:length(chr_order)){
  a<-Gene.Label[Gene.Label$Chromosome==chr_order[i],]%>%arrange(po)
  b<-rbind(b,a)
}
Gene.Label<-b
colnames(Gene.Label)[5]<-"chr"

rm(a,b)

map1<-Gene.Label[,c(5,10,2)]%>%distinct(gene,.keep_all = T)
colnames(map1)<-c("chr","po","Gene")
map2<-Gene.Label[,c(5,10,1)]%>%distinct(Probe,.keep_all = T)
colnames(map2)<-c("chr","po","Probe");
colnames(probe_exp_meth_cr)
methy<-probe_exp_meth_cr[,c(13:14,7,1)]%>%distinct(Probe,.keep_all = T) %>% 
  filter(Probe%in%probeset1$Probe) %>% 
  mutate(value=ifelse(logFC_m_tcga>0,1,-1),
         Chromosome=substring(Chromosome,4,)) %>% 
  dplyr::select(-c("logFC_m_tcga","Probe"))
methy<-map2 %>% 
  inner_join(probe_exp_meth_cr[,c(1,7)],by="Probe")%>% 
  mutate(value=ifelse(logFC_m_tcga>0,1,-1),
         color=ifelse(logFC_m_tcga>0,"red","blue"))%>% 
  dplyr::select(-"logFC_m_tcga") 
table(probeset1$logFC_m_tcga>0)
table(probe_exp_meth_cr$logFC_m_tcga>0)



colnames(methy)[3]<-"NAME"
colnames(methy)<-c("chr","start","end","value")
exp<-map2 %>% 
  inner_join(probe_exp_meth_cr[,c(1,8)],by="Probe")%>% 
  mutate(value=ifelse(logFC_g_tcga>0,1,-1),
         color=ifelse(logFC_g_tcga>0,"red","blue"))%>% 
  dplyr::select(-"logFC_g_tcga")  
colnames(exp)[3]<-"NAME"

library(OmicCircos)
pdf("./Methy_Diff/Circos.pdf",width=8,height=8)
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
circos(R=180, type="chr", cir="hg19", print.chr.lab=T, W=4);
circos(R=170, cir="hg19", W=20, mapping=map1, type="label", 
       side="in", col="black",cex=0.4);
circos(R=200, cir="hg19", W=20, mapping=map2, type="label", 
       side="out",  cex=0.6);

dev.off()

pdf("./Methy_Diff/Circos.pdf",width=8,height=8)
par(mar=c(0, 0, 0, 0));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="");
circos(R=200, type="chr", cir="hg19", print.chr.lab=T, W=4);
circos(R=220, cir="hg19", W=20, mapping=map1, type="label", 
       side="out", col="black",cex=0.6);
colors<-Gene.Label$cr;
circos(R=180, cir="hg19", W=20, mapping=map2, type="label", 
       side="in",cex=0.6,col=colors);

dev.off()
## ImPS related gene----
rm(list=ls())
load("/Volumes/T7/Methylation/Download/06Imps_tcga_sp.rda")
load("/Volumes/T7/Methylation/GEO/GSE66863_GSE66836_Data/06Imps_gse66836_sp.rda")
load("/Volumes/T7/Methylation/GEO/GSE60644_GSE60645_Data/06Imps_gse56044_sp.rda")
load("/Volumes/T7/Methylation/00Gene_exp_methy_match.rda")
probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe,
                probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
gene<-unique(c(probe_exp_methy_d_cis$gene,probe_exp_methy_u_cis$gene,
               probe_exp_methy_d_trans$gene,probe_exp_methy_u_trans$gene))
probe_exp_methy<-rbind(probe_exp_methy_d_cis,
                       probe_exp_methy_d_trans,
                       probe_exp_methy_u_cis,
                       probe_exp_methy_u_trans)

imps_up<-Reduce(intersect,list(gse56044_imps[gse56044_imps$Imps_Correlation>0,"Probe"],
                               gse66836_imps[gse66836_imps$Imps_Correlation>0,"Probe"],
                               tcga_imps[tcga_imps$Imps_Correlation>0,"Probe"]))
imps_down<-Reduce(intersect,list(gse56044_imps[gse56044_imps$Imps_Correlation<0,"Probe"],
                                 gse66836_imps[gse66836_imps$Imps_Correlation<0,"Probe"],
                                 tcga_imps[tcga_imps$Imps_Correlation<0,"Probe"]))
gse56044_imps<-gse56044_imps%>%column_to_rownames("Probe")
rownames(tcga_imps)<-tcga_imps$Probe
rownames(gse66836_imps)<-gse66836_imps$Probe

imps_up<-cbind(Probe=imps_up,
               gse56044_imps=gse56044_imps[imps_up,],
               tcga_imps=tcga_imps[imps_up,2],
               gse66836_imps=gse66836_imps[imps_up,2])
imps_down<-cbind(Probe=imps_down,
                 gse56044_imps=gse56044_imps[imps_down,],
                 tcga_imps=tcga_imps[imps_down,2],
                 gse66836_imps=gse66836_imps[imps_down,2])
imps<-rbind(imps_up,imps_down)%>%as.data.frame()%>%
  inner_join(probe_exp_methy,by="Probe")




load("/Volumes/T7/Methylation/Download/06Imps_exp_tcga_sp.rda")
load("/Volumes/T7/Methylation/GEO/GSE66863_GSE66836_Data/06Imps_exp_gse66836_sp.rda")
load("/Volumes/T7/Methylation/GEO/GSE60644_GSE60645_Data/06Imps_exp_gse56044_sp.rda")
imps_up_g<-Reduce(intersect,list(gse56044_exp_imps[gse56044_exp_imps$Imps_Correlation>0,"Gene"],
                                 gse66836_exp_imps[gse66836_exp_imps$Imps_Correlation>0,"Gene"],
                                 tcga_exp_imps[tcga_exp_imps$Imps_Correlation>0,"Gene"]))
imps_down_g<-Reduce(intersect,list(gse56044_exp_imps[gse56044_exp_imps$Imps_Correlation<0,"Gene"],
                                   gse66836_exp_imps[gse66836_exp_imps$Imps_Correlation<0,"Gene"],
                                   tcga_exp_imps[tcga_exp_imps$Imps_Correlation<0,"Gene"]))

gse56044_exp_imps<-gse56044_exp_imps%>%column_to_rownames("Gene")
rownames(tcga_exp_imps)<-tcga_exp_imps$Gene
rownames(gse66836_exp_imps)<-gse66836_exp_imps$Gene

imps_up_g<-cbind(gene=imps_up_g,
                 gse56044_imps_g=gse56044_exp_imps[imps_up_g,],
                 tcga_imps_g=tcga_exp_imps[imps_up_g,2],
                 gse66836_imps_g=gse66836_exp_imps[imps_up_g,2])
imps_down_g<-cbind(gene=imps_down_g,
                   gse56044_imps_g=gse56044_exp_imps[imps_down_g,],
                   tcga_imps_g=tcga_exp_imps[imps_down_g,2],
                   gse66836_imps_g=gse66836_exp_imps[imps_down_g,2])
imps_g<-rbind(imps_up_g,imps_down_g)%>%as.data.frame()
imps_Pheno<-imps[imps$gene%in%imps_g$gene,]%>%
  inner_join(imps_g,by="gene") %>%
  dplyr::select(c(1,5,12:18,2:4,19:21,6:11))
load("./00Probe_filtercr_sp.rda")
x<-imps_Pheno[imps_Pheno$Probe%in%intersect(probe_exp_meth_cr$Probe,imps_Pheno$Probe),]
colnames(imps_Pheno)
save(imps_Pheno,file="./01IMps_Pheno_Probe_sp.rda")
load("./01IMps_Pheno_Probe_sp.rda")
write.csv(imps_Pheno,file="./TableS3IMPS_pheno.csv")

probeset1<-intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe)
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%probeset1,"gene"])


## mIMg----
rm(list=ls())
load("./01IMps_Pheno_Probe_sp.rda")
load("./00Probe_filtercr_sp.rda")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"Probe"]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"gene"])
result<-probe_exp_meth_cr[probe_exp_meth_cr$Probe%in%probeset1,c(1,6,13,14,15,16,17,18,3,2,7,8,4,9,10,5,11,12)]
result$methylation_status<-ifelse(result$logFC_m_tcga>0,"Hyper","Hypo")
result$expression_status<-ifelse(result$logFC_g_tcga>0,"Up-regulated","Down-regulated")
result<-result[,c(1:9,19:20,10:18)]
colnames(probeset1)
result<-result %>% 
  inner_join(imps_Pheno[,c(1,11,14,10,13,12,15)],by="Probe")
result$`methylation&immune`<-ifelse(result$tcga_imps>0,"Cis","Trans")
result$`expression&immune`<-ifelse(result$tcga_imps_g>0,"Cis","Trans")
result<-result[,c(1:20,27,21,23,25,28,22,24,26)]
write.csv(result,file="Table_S6mIMg.csv")
save(result,file="./00mIMg.rda")
## Figure 4 D-E-------
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
gene.df <- bitr(geneset1,
                fromType = "SYMBOL",
                toType = c( "ENTREZID"),
                OrgDb = org.Hs.eg.db)#
m_tg <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
em<-enricher(geneList$ENTREZID,TERM2GENE = m_tg)
em2 <- GSEA(gene.df$ENTREZID, TERM2GENE = m_tg)
head(em)
head(m_t2g)
geneList<-imps_Pheno[imps_Pheno$gene%in%geneset1,c("gene","logFC_g_tcga")] %>%
  distinct(gene,.keep_all = T) %>% 
  rename("SYMBOL"="gene") %>% 
  inner_join(gene.df[,c("SYMBOL","ENTREZID")],by="SYMBOL") %>% 
  arrange(desc(logFC_g_tcga))
ego<- enrichGO(gene = geneList$ENTREZID,ont   = "BP",
               OrgDb = org.Hs.eg.db, 
               readable = TRUE)

edox <- setReadable(ego ,'org.Hs.eg.db', 'ENTREZID')
pdf("./Phenotype/mIMg.pdf",height=6,width=10)
cnetplot(edox, circular = TRUE, colorEdge = TRUE,cex_label_category=0)
dev.off()
ego<- enrichGO(gene = geneList$ENTREZID,ont   = "MF",
               OrgDb = org.Hs.eg.db, 
               readable = TRUE)
edox <- setReadable(ego,'org.Hs.eg.db', 'ENTREZID')
pdf("./Phenotype/mIMg2.pdf",height=6,width=8)
cnetplot(edox,foldChange=geneList$logFC_g_tcga)
dev.off()


### Figure 4 B Cytoscape link------
library(networkD3)
Links<-data.table::fread("./network/string_interactions.tsv")
probeset1<-imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"Probe"]
geneset1<-unique(imps_Pheno[imps_Pheno$Probe%in%intersect(imps_Pheno$Probe,probe_exp_meth_cr$Probe),"gene"])
colnames(probe_exp_meth_cr)
result<-probe_exp_meth_cr[probe_exp_meth_cr$Probe%in%probeset1,c(1,2,6,3,7,8)]
result$methylation_status<-ifelse(result$logFC_m_tcga>0,"Hyper","Hypo")
result$expression_status<-ifelse(result$logFC_g_tcga>0,"Up-regulated","Down-regulated")
result<-result %>% 
  inner_join(imps_Pheno[,c(1,11,14)],by="Probe")
result$`methylation&immune`<-ifelse(result$tcga_imps>0,"Cis","Trans")
result$`expression&immune`<-ifelse(result$tcga_imps_g>0,"Cis","Trans")
Mislinks<-result[,c("Probe","gene","Estimate_tcga","cr")]
Mislinks$type<-"probe_gene"
colnames(Mislinks)[1:3]<-c("source","target","value")
colnames(Links)
Mislinks2<-Links[,c(1,2,13)]
colnames(Mislinks2)<-c("source","target","value")
Mislinks2$cr<-NA
Mislinks2$type<-"gene_gene"
Mislinks$value<-abs(Mislinks$value)
Mislinks<-rbind(Mislinks,Mislinks2)
write.table(Mislinks,"./network/mislinks.txt",sep="\t",quote=F,col.names = T,row.names = F)
Node<-data.frame(Node=c(result$Probe,unique(result$gene)),
                 type=c(rep("Probe",length(result$Probe)),rep("gene",length(unique(result$gene)))),
                 imps=c(result$tcga_imps,unique(result$tcga_imps_g)))
write.table(Node,"./network/node.txt",sep="\t",quote=F,col.names = T,row.names = F)


### 