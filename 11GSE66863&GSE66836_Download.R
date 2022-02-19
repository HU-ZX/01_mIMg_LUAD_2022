# GSE66863 GPL14550
rm(list=ls())
library(tidyverse)
library("GEOquery")
files<-list.files("GSE66863_GSE66836_Data/GSE66863", full.names = T)
target<-getGEO(GEO="GSE66863",file=files[12],getGPL=F)
pda<-pData(target)
exp<-exprs(target)
colnames(pda)
pda[1,]
### Annotation_GPL----
path<-"./GSE66863/GPL14550_reannotation"
dir.create(path)
library(limma)
anno<-data.table::fread(files[4])
anno_seqs <- anno %>% 
  select("ID", "SEQUENCE") %>% 
  mutate(seqs = paste0(">", ID, "\n", SEQUENCE))
anno_seqs<-anno_seqs[12:nrow(anno_seqs),]
write_lines(anno_seqs$seqs, path = paste0(path,"/GPL14550.fa"))
getwd()
library(Rsubread)
align(index = "/Volumes/T7/Index_Gencode_v35/Gencode_v35",
      readfile1 = paste0(path,"/GPL14550.fa"),
      output_file = paste0(path,"/GPL14550_gencode.BAM"),
      nBestLocations = 16,
      keepReadOrder = T)
library(Rsamtools)
bam_res <- scanBam(file = paste0(path,"/GPL14550_gencode.BAM")) %>% 
  .[[1]] %>% 
  .[1:11] %>% 
  bind_cols() %>% 
  separate(col = rname, into = paste0("id", 1:8),sep = "\\|",remove = T) %>% 
  add_count(cigar, name = "cigar_n") %>% 
  filter(cigar == "60M") 

p2s_rsubread <- bam_res %>% 
  select(probe_id = qname, symbol = id6) %>% 
  distinct() %>% 
  add_count(probe_id) %>% 
  filter(n == 1) %>% 
  select(1:2)
write_tsv(p2s_rsubread, 
          path = paste0(path,"/GPL14550_reannotation.tsv"))
### pDATA----
# The gene expression values were log 2 transformed 
# and normalized between arrays by using the 75th percentile 
# method in Genespring GX analysis Software v.12.1 (Agilent technology).
pda<-pda[,c(1:2,10:17)]
pda[1,]
colnames(pda)[3:10]<-c("Age","Smoking_history","Gender",
                       "Pathology","Stage","Tp53_status",
                       "EGFR_status","KRAS_status")
table(pda$Pathology)
pda<-pda%>%
  mutate(Age=as.numeric(sapply(str_split(Age,":"),"[",2)),
         Smoking_history=substring(sapply(str_split(Smoking_history,":"),"[",2),2,),
         Gender=substring(sapply(str_split(Gender,":"),"[",2),2,),
         Stage=substring(sapply(str_split(Stage,":"),"[",2),2,),
         Tp53_status=substring(sapply(str_split(Tp53_status,":"),"[",2),2,),
         EGFR_status=substring(sapply(str_split(EGFR_status,":"),"[",2),2,),
         KRAS_status=substring(sapply(str_split(KRAS_status,":"),"[",2),2,))%>%
  dplyr::select(-"Pathology")
files_sup<-data.frame(Files=list.files("./GSE66863/GSE66863_RAW", full.names = F))%>%
  mutate(geo_accession=str_extract(Files,"GSM\\d+"))
pda<-pda%>%inner_join(files_sup,by="geo_accession")

### Exprs___low quality----
exp<-exprs(target)%>%as.data.frame()%>%
  rownames_to_column("probe_id")%>%
  inner_join(p2s_rsubread,by="probe_id")%>%
  dplyr::select(-"probe_id")%>%
  aggregate(.~symbol,data=.,max)%>%
  mutate(gene_symbol=symbol)%>%
  dplyr::select(-"symbol")%>%
  dplyr::select("gene_symbol",everything())
exp2<-exp%>%column_to_rownames("gene_symbol")%>%as.matrix()
boxplot(log2(exp2),
        ylab = expression(log[2](intensity)),
        las = 2,
        outline = FALSE)
## read.maimages----
pre<-data.table::fread(paste0("./GSE66863/GSE66863_RAW/",pda$Files[1]),skip=9)
?fread()
GSE66863_raw <- read.maimages(files = pda$Files,
                              source = "agilent",
                              path = "./GSE66863/GSE66863_RAW",
                              names = pda$geo_accession,
                              other.columns = "gIsWellAboveBG",
                              green.only = T)

colnames(pre)
anno<-pre[,13:14]
?read.maimages()
GSE66863_raw$targets<-pda
table(GSE66863_raw$genes$ControlType)
head(GSE66863_raw$E)[, 1:5]
dim(GSE66863_raw)
## Quality Control----
boxplot(log2(GSE66863_raw$E),
        ylab = expression(log[2](intensity)),
        las = 2,,
        outline = FALSE)
library(arrayQualityMetrics)
rownames(pda) <- pda$geo_accession
identical(colnames(GSE66863_raw$E),pda$geo_accession)
eset<-GSE66863_raw$E
eset <- ExpressionSet(assayData = GSE66863_raw$E, 
                      phenoData = AnnotatedDataFrame(data = pda))

arrayQualityMetrics(eset, 
                    outdir = "./GSE66863/GSE66863_QC_raw", 
                    do.logtransform = TRUE,
                    force = TRUE)
dev.off()


dir.create("./GSE66863/GSE66863_QC_raw/image_plots/")
for(i in ncol(GSE66863_raw$E)){
  tiff(filename = paste0("./GSE66863/GSE66863_QC_raw/image_plots/",
                         pda$geo_accession[i], "_", i,".tif"),
       height = 10, 
       width = 10, 
       units = "in",
       res = 300)
  image(log2(GSE66863_raw$E[, i]),GSE66863_raw$printer,main = pda$geo_accession[i])
  dev.off()
}
### Normalization----
anno<-data.table::fread("./GSE66863/GPL14550_reannotation/GPL14550_reannotation.tsv")
bgc <-  backgroundCorrect(RG = GSE66863_raw, 
                          method = "normexp",
                          offset = 50,
                          normexp.method = "mle")

norm <- normalizeBetweenArrays(bgc, method = "quantile")
nrow(expr_norm)
sum(is.na(expr_norm))
rownames(expr_norm)<-norm$genes$GeneName
norm$genes$GeneName[duplicated(norm$genes$GeneName)]
norm$genes$GeneName[str_detect(norm$genes$GeneName,"HLA")]
expr_norm <- norm$E%>%as.data.frame()%>%
  mutate(Gene=norm$genes$GeneName)%>%
  aggregate(.~Gene,data=.,max)

# expr_norm2 <- norm$E%>%as.data.frame()%>%
#   mutate(probe_id=norm$genes$ProbeName)%>%
#   inner_join(anno,by="probe_id")%>%
#   mutate(Gene=symbol)%>%
#   dplyr::select(-c("probe_id","symbol"))%>%
#   aggregate(.~Gene,data=.,max)
rt<-expr_norm%>%column_to_rownames("Gene")%>%
  as.matrix()
boxplot(norm$E,
        ylab = expression(log[2](intensity)),
        las = 2,,
        outline = FALSE)
boxplot(rt,
        ylab = expression(log[2](intensity)),
        las = 2,,
        outline = FALSE)


save(norm,pda,expr_norm,GSE66863_raw,file="./01Importexpression.Rda")



