#### 01 Download IDAT Files----
{
  projects <- TCGAbiolinks:::getGDCprojects()$project_id
  projects <- projects[grepl('^TCGA',projects,perl=T)]
  match.file.cases.all <- NULL
  query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Raw microarray data",
                  data.type = "Raw intensities", 
                  experimental.strategy = "Methylation array", 
                  legacy = TRUE,
                  file.type = ".idat",
                  platform = "Illumina Human Methylation 450")
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- "TCGA-LUAD"
  match.file.cases.all <- match.file.cases
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
        error = function(e) GDCdownload(query, method = "client"))
  }
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
readr::write_csv(match.file.cases.all, path =  "idat_filename_case.csv")

#### 02 Prepare Target File ----
rm(list=ls())
getwd()
files <- list.files("/Volumes/T7/Methylation/Download/GDCdata/TCGA-LUAD/legacy/Raw_microarray_data/Raw_intensities", full.names = T)
n <- length(files)
head(files)
new_dir <- paste(getwd(),"/merge_dir",sep = "")
dir.create(new_dir)
for(i in 1:n){
  aaa =list.files(files[i],pattern = ".idat",full.names=T)
  file.copy(from = aaa,to = new_dir)
}
# x<-list.files(files,pattern = ".idat")
# x1<-list.files(files,pattern = ".idat",full.names=T)
# table(is.na(x))
# y<-list.files(new_dir,pattern = ".idat")
# setdiff(x,y)
# library(tidyverse)
# x1[str_detect(x1,setdiff(x,y))]
# length(list.files(new_dir))

library("minfi")
Filename<-list.files(new_dir,pattern = ".idat")
str_locate_all(Filename[1],"\\_")
Basename<-substring(Filename,1,17)
Array<-sapply(strsplit(Basename,"\\_"),"[",2)
Slide<-sapply(strsplit(Basename,"\\_"),"[",1)
target_sheet<-data.frame(Basename=Basename,
                         Array=Array,
                         Slide=Slide)%>%
  distinct(Basename,.keep_all = T)
# targets <- read.metharray.sheet(getwd())
targets <- data.table::fread("./idat_filename_case.txt")%>%
  mutate(Basename=substring(file_name,1,17))%>%
  distinct(Basename,.keep_all = T)%>%
  dplyr::select(-c("file_name"))%>%
  inner_join(target_sheet,by="Basename")%>%
  mutate(Group=substring(cases,14,16),
         Sample=substring(cases,1,12),
         Basename=paste0(new_dir,"/",Basename))%>%
  # dplyr::filter(!(Group%in%c("11A","02A")))%>%
  dplyr::filter(Group=="01A")%>%
  distinct(Sample,.keep_all = T)
  

table(targets$Group) 
  
dup_sample <- targets %>%
  add_count(Sample)%>%
  filter(!(n == 1))%>%
  arrange(Sample)%>%
  filter(!(Group=="11A"))
Time<-data.table::fread("/Volumes/T7/Pan_cancer_Clinic/Pancancer_survivalsup.txt")
Clinic<-data.table::fread("/Volumes/T7/Methylation/GDC_clinical_cart/clinical.tsv")%>%
  distinct(case_submitter_id,.keep_all = T)%>%
  dplyr::select(c(1,2,4,10,12,16,25:28,48,108,110,116))
vacant<-Clinic$days_to_death[1]
Clinic$days_to_death[str_detect(Clinic$days_to_death,vacant)]<-NA
Clinic$age_at_index[str_detect(Clinic$age_at_index,vacant)]<-NA
Clinic$days_to_last_follow_up[str_detect(Clinic$days_to_last_follow_up,vacant)]<-NA
Clinic$days_to_death[is.na(Clinic$days_to_death)]<-Clinic$days_to_last_follow_up[is.na(Clinic$days_to_death)]
Clinic<-Clinic[,-11]
colnames(Clinic)<-c("case_id","Sample","Age","OS.Time","Gender","OS_Vital_Status",
                    "M_Stage","N_Stage","Stage","T_Stage","Primary_Diagnosis","Prior_Malignancy","Site_of_Biopsy")  
Time<-Time[Time$type=="LUAD",]%>%
  dplyr::select(c(2,4,5))
colnames(Time)<-c("Sample","PFI","PFI.time")
Clinic<-Clinic%>%inner_join(Time,by="Sample")
table(Clinic$OS_Vital_Status)
Clinic$OS_Vital_Status<-ifelse(Clinic$OS_Vital_Status=="Alive","0","1")
Clinic$OS_Vital_Status<-as.numeric(Clinic$OS_Vital_Status)
Clinic$OS.Time<-as.numeric(Clinic$OS.Time)
table(Clinic$Stage)
Clinic$Stage[str_detect(Clinic$Stage,"III")]<-"Stage III"
Clinic$Stage[Clinic$Stage%in%c("Stage II","Stage IIA","Stage IIB")]<-"Stage II"
Clinic$Stage[Clinic$Stage%in%c("Stage I","Stage IA","Stage IB")]<-"Stage I"
Clinic$Stage[str_detect(Clinic$Stage,vacant)]<-NA
table(Clinic$M_Stage)
Clinic$M_Stage[str_detect(Clinic$M_Stage,vacant)]<-"MX"
Clinic$M_Stage[str_detect(Clinic$M_Stage,"M1")]<-"M1"
table(Clinic$T_Stage)
Clinic$T_Stage[str_detect(Clinic$T_Stage,"T1")]<-"T1"
Clinic$T_Stage[str_detect(Clinic$T_Stage,"T2")]<-"T2"
table(Clinic$N_Stage)
Clinic$N_Stage[str_detect(Clinic$N_Stage,vacant)]<-"NX"
table(Clinic$Prior_Malignancy)
table(Clinic$Site_of_Biopsy)
table(Clinic$Gender)
table(Clinic$Primary_Diagnosis)
a<-intersect(Clinic$Sample,targets$Sample)
Targets<-targets%>%
  inner_join(Clinic,by="Sample")
table(duplicated(targets$Sample))
write.table(Targets,"Targets.txt",quote=F,sep="\t",col.names = T,row.names = F)

#### 03 Download Expression-----
rm(list=ls())
table<-TCGAbiolinks::getSampleFilesSummary("TCGA-LUAD")
query.exp <- GDCquery(project = "TCGA-LUAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts")
??GDCquery()
GDCdownload(query.exp)

data <- GDCprepare(query.exp, save = TRUE, 
                   save.filename = "luad_expression.rda",
                   remove.files.prepared = FALSE)
rowRanges(data)
Matrix <- assay(data)

json_file <- './TCGALUAD_metadata.cart.2021-02-16.json'
metadata <- jsonlite::read_json(path = json_file, simplifyVector = T)
metadata <- tibble::tibble(
  file_name = metadata$file_name,
  md5sum = metadata$md5sum,
  TCGA_id_full = bind_rows(metadata$associated_entities)$entity_submitter_id,
  TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
  patient_id = stringr::str_sub(TCGA_id, 1, 12),
  tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
  tissue_type = sapply(tissue_type_id, function(x) {
    switch(x,
           "01" = "Primary Solid Tumor",
           "02" = "Recurrent Solid Tumor",
           "03" = "Primary Blood Derived Cancer - Peripheral Blood",
           "05" = "Additional - New Primary",
           "06" = "Metastatic",
           "07" = "Additional Metastatic",
           "11" = "Solid Tissue Normal")}),   
  group = ifelse(tissue_type_id == "11", "Normal", "Tumor"))
Path = "/Volumes/T7/Methylation/Download/GDCdata/TCGA-LUAD/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"
files = list.files(Path,full.names = T)
n = length(files)
new_dir = paste(getwd(),"/count_dir",sep = "")
dir.create(paste(getwd(),"/count_dir",sep = ""))
for(i in 1:n){
  aaa = list.files(files[i],pattern = "htseq.counts.gz",full.names = T)
  file.copy(from = aaa,to = new_dir)
}

counts_files <- list.files(path = new_dir, 
                           pattern = "htseq.counts.gz", 
                           full.names = T, 
                           recursive = T)
test <- readr::read_tsv(counts_files[1], col_names = F)
Targets<-data.table::fread("./Targets.txt")
require("tools")
all(tools::md5sum(counts_files) %in% metadata$md5sum)
purrr::set_names("gene_id", basename(counts_files[1]))
counts_df <- counts_files %>% 
  lapply(function(x) {
    tmp <- read_tsv(x, col_names = F) %>% 
      purrr::set_names("gene_id", basename(x))
    cat(which(counts_files == x), "of", length(counts_files), "\n")
    return(tmp)
  }) %>%
  reduce(function(x, y) full_join(x, y, by = "gene_id")) %>% 
  dplyr::select(gene_id, metadata$file_name) %>% 
  set_names("gene_id", metadata$TCGA_id_full) %>% 
  dplyr::slice(1:(nrow(.)-5))
getwd()
load("/Volumes/T7/Index_Gencode_v35/anno.Rdata")
dim(anno)#60612

metadata<-metadata%>%
  filter(group=="Tumor")%>%
  filter(substring(TCGA_id,14,16)=="01A")%>%
  arrange(TCGA_id_full)%>%
  distinct(patient_id,.keep_all = T)%>%
  mutate(Sample=patient_id)
table(metadata$tissue_type)
Targets1<-Targets%>%
  inner_join(metadata,by="Sample")
table(substr(colnames(counts_df ),14,15)=="01") 
count<-counts_df%>%
  mutate(gene_id=substring(counts_df$gene_id,1,15))%>%
  column_to_rownames("gene_id")%>%
  dplyr::select(Targets1$TCGA_id_full)
overlap<-anno$ENSG%in%rownames(count)
anno0<-anno[overlap,]
overlap<-rownames(count)%in%anno$ENSG
count0<-count[overlap,]
R<-apply(count0,2,sum)
rownames(anno0)<-anno0$ENSG
flg<-anno0[rownames(count0),'Length']

TPM<-count0
for(i in 1:ncol(TPM)){
  t<-TPM[,i]/flg
  T<-sum(t)
  TPM[,i]<-t*10^6/T
}
colnames(TPM)<-substring(colnames(TPM),1,12)
identical(colnames(TPM),Targets1$patient_id)
TPM_sym<-TPM%>%
  mutate(ENSG=rownames(TPM))%>%
  inner_join(anno[,4:5],by="ENSG")%>%
  distinct(Symbol,.keep_all = T)%>%
  column_to_rownames("Symbol")%>%
  dplyr::select(-"ENSG")
count<-count0
save(counts_df,count0,file="00.OriExp.Rda")
Targets<-Targets%>%
  filter(Sample%in%Targets1$Sample)
Targets<-Targets[,-"case_id"]
save(TPM_sym,TPM,Targets,file="02.TPM&Targets.Rda")