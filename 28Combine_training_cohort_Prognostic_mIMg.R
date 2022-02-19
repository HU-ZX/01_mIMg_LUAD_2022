# Get prognostic mIMg probeset4 geneset4----
load("/Users/huzixin/Data_Project/Methylation/00mIMg.rda")
mIMg<-result
load("/Users/huzixin/Data_Project/Methylation/Validation/Validation_Result/gpl_multicox.rda")
load("/Users/huzixin/Data_Project/Methylation/Download/tcga_Survival_Results/tcga_multicox1_probe.rda")
load("/Users/huzixin/Data_Project/Methylation/Download/tcga_Survival_Results/tcga_multicox1.rda")

geneset<-mIMg %>% 
  filter(Probe%in%result1_tcga_probe$Probe[result1_tcga_probe$P.val<0.05]) %>%
  filter(gene%in%result1_tcga$gene[result1_tcga$P.val<0.05]) %>% 
  dplyr::select("gene") %>% 
  distinct(gene)

probeset<-mIMg %>% 
  filter(Probe%in%result1_tcga_probe$Probe[result1_tcga_probe$P.val<0.05]) %>%
  filter(gene%in%result1_tcga$gene[result1_tcga$P.val<0.05])  %>% 
  arrange(gene) %>% 
  dplyr::select(-c(9:21)) %>% 
  inner_join(result1_tcga_probe[,-6],by="Probe")


# Multicox Table1 ----

result<-probeset[,c("Probe","gene")] %>% 
  distinct(gene,.keep_all = T) %>% 
  inner_join(result_gpl570,by="gene") %>% 
  inner_join(result_gpl6884,by="gene") %>% 
  inner_join(result1_tcga,by="gene") %>% 
  filter(P.val<0.1&P.val.x<0.1&P.val.y<0.1)%>% 
  filter((HR>1&HR.x>1&HR.y>1)|(HR<1&HR.x<1&HR.y<1)) %>% 
  arrange(gene)
#Table1
write.table(result,"./02all_multicox.txt",quote=F,sep="\t",col.names=T,row.names=F)                                


colnames(imps_Pheno)
result_probe<-result1_tcga_probe %>% 
  inner_join(mIMg[,c(1,3:8)],by="Probe") %>% 
  #filter(Probe%in%probe_exp_meth_cr$Probe) %>% 
  filter(Probe%in%probeset$Probe) %>% 
  #filter(gene%in%result$gene) %>% 
  #inner_join(result,by="gene") %>% 
  arrange(gene)

write.table(result_probe,"./02all_multicox_probe.txt",quote=F,sep="\t",col.names=T,row.names=F)                                


colnames(result_gpl570)<-paste0(c("gene","HR","group","95%CI","P.val"),"_gpl570")
colnames(result_gpl570)[1]<-"gene"
colnames(result_gpl6884)<-paste0(c("gene","HR","group","95%CI","P.val"),"_gpl6884")
colnames(result_gpl6884)[1]<-"gene"
vali_result<-geneset4 %>% 
  inner_join(result1_tcga,by="gene") %>% 
  inner_join(result_gpl570,by="gene") %>%
  left_join(result_gpl6884,by="gene") 
write.table(vali_result,"./02all_multicoxvali.txt",quote=F,sep="\t",col.names=T,row.names=F)                                

# Unicox----                                
load("/Users/huzixin/Data_Project/Methylation/Validation/Validation_Result/gpl_unicox.rda")
load("/Users/huzixin/Data_Project/Methylation/Download/tcga_Survival_Results/tcga_unicox1.rda")                               
colnames(Unicox_result_gpl570)[2:4]<-paste0(colnames(Unicox_result_gpl570)[2:4],"_gpl570")
colnames(Unicox_result_gpl6884)[2:4]<-paste0(colnames(Unicox_result_gpl6884)[2:4],"_gpl6884")

Uni_vali_result<-geneset4 %>% 
  inner_join(Unicox_result1_tcga,by="gene") %>% 
  inner_join(Unicox_result_gpl570,by="gene") %>%
  left_join(Unicox_result_gpl6884,by="gene") 
write.table(Uni_vali_result,"./02all_unicoxvali.txt",quote=F,sep="\t",col.names=T,row.names=F)                                
