rm(list=ls())
library(tidyverse)
load("/Users/huzixin/Data_Project//Methylation/00Gene_exp_methy_match.rda")
load("./01Importexpression.Rda")
load("./02Methylation_byChAMP.Rda")
GSE66836_Methy<-data.table::fread("./GSE66836_myDMP.txt")
cis_probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe))
trans_probe<-unique(c(probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
rownames(pD)<-pD$geo_accession_m
pD<-pD[colnames(beta_norm),]
identical(colnames(beta_norm),rownames(pD))
colnames(beta_norm)<-pD$title
rownames(pD)<-pD$geo_accession_e
expr_norm<-expr_norm%>%
  column_to_rownames("Gene")%>%
  dplyr::select(pD$geo_accession_e)%>%
  as.matrix()
identical(colnames(expr_norm),rownames(pD))
colnames(expr_norm)<-pD$title
gse66836_exp<-expr_norm
gse66836_beta<-beta_norm
identical(colnames(gse66836_exp),colnames(gse66836_beta))
rm(beta_norm,expr_norm,GSE66863_raw,Myload,norm,pda)
save(gse66836_beta,gse66836_exp,file="./04exp66836_methy66863_data")
probe<-cis_probe
probe_fl<-c()
esm<-c()
p.val<-c()
for(i in 1:length(probe)){
  a<-gse66836_beta[probe[i],]
  gene<-GSE66836_Methy[GSE66836_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(gse66836_exp)){
    b<-gse66836_exp[gene,]%>%as.numeric()
  }
  if(!anyNA(b)){
    test<-cor.test(a,b,method="spearman")
    if(!anyNA(c(test$p.value,test$estimate))){
      if(test$p.value<0.05 & abs(test$estimate)>0.5){
        probe_fl<-c(probe_fl,probe[i]) 
        esm<-c(esm,test$estimate)
        p.val<-c(p.val,test$p.value)
      } 
    }
  }
}
probe_fl_66836_cis<-probe_fl
esm_66836_cis<-esm
pval_66836_cis<-p.val

probe_fl_df_66836_cis<-data.frame(Probe=probe_fl_66836_cis,
                                  Estimate=esm_66836_cis,
                                  P.val=pval_66836_cis)
probe_fl_df_66836_cis$cr<-rep("Cis",nrow(probe_fl_df_66836_cis))
probe<-trans_probe
probe_fl<-c()
esm<-c()
p.val<-c()

for(i in 1:length(probe)){
  a<-gse66836_beta[probe[i],]
  gene<-GSE66836_Methy[GSE66836_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(gse66836_exp)){
    b<-gse66836_exp[gene,]%>%as.numeric()
  }
  if(!anyNA(b)){
    test<-cor.test(a,b,method="spearman")
    if(!anyNA(c(test$p.value,test$estimate))){
      if(test$p.value<0.05 & abs(test$estimate)>0.5){
        probe_fl<-c(probe_fl,probe[i]) 
        esm<-c(esm,test$estimate)
        p.val<-c(p.val,test$p.value)
      } 
    }
  }
}

probe_fl_66836_trans<-probe_fl
esm_66836_trans<-esm
pval_66836_trans<-p.val

probe_fl_df_66836_trans<-data.frame(Probe=probe_fl_66836_trans,Estimate=esm_66836_trans,P.val=pval_66836_trans)
probe_fl_df_66836_trans$cr<-rep("Trans",nrow(probe_fl_df_66836_trans))
probe_fl_66836<-rbind(probe_fl_df_66836_cis,probe_fl_df_66836_trans)
colnames(probe_fl_66836)[2:3]<-paste0(colnames(probe_fl_66836)[2:3],"_gse66836")
save(probe_fl_66836,file="04Co_exp_methy_probe_fl_sp.rda")
# combine logFC-----
rm(list=ls())
load("./05Diffgene_gse66863.rda")
table(Df_gse66863$logFC<0)
table(Df_gse66863$logFC>0)
GSE66836_Methy<-data.table::fread("./GSE66836_myDMP.txt")
table(GSE66836_Methy$logFC<0)
table(GSE66836_Methy$logFC>0)
Methy<-GSE66836_Methy%>%select("Probe","logFC","gene")
Expr<-Df_gse66863%>%select("Symbol","logFC")
colnames(Expr)[c(2,1)]<-c("logFC_g","gene")
colnames(Methy)[2]<-"logFC_m"

gse66836_Expr_Methy<-Expr%>%
  inner_join(Methy,by="gene")%>%
  dplyr::select("Probe","gene","logFC_g","logFC_m")
cis_downgene<-gse66836_Expr_Methy[gse66836_Expr_Methy$logFC_g<0&gse66836_Expr_Methy$logFC_m<0,]
cis_upgene<-gse66836_Expr_Methy[gse66836_Expr_Methy$logFC_g>0&gse66836_Expr_Methy$logFC_m>0,]
trans_downgene<-gse66836_Expr_Methy[gse66836_Expr_Methy$logFC_g<0&gse66836_Expr_Methy$logFC_m>0,]
trans_upgene<-gse66836_Expr_Methy[gse66836_Expr_Methy$logFC_g>0&gse66836_Expr_Methy$logFC_m<0,]


save(gse66836_Expr_Methy,file="05exp_methy_logFC.rda")

# Figure 3 G Volcano----
library(tidyverse)
rm(list=ls())
load("./05Diffgene_gse66863.rda")
gse66863_Methy<-data.table::fread("./GSE66836_myDMP.txt")
load("/Users/huzixin/Data_Project/Methylation/00Gene_exp_methy_match.rda")
colnames(gse66863_Methy)
colnames(Df_gse66863) 
methy<-gse66863_Methy[,c(1,2,6,15)]
exp_g<-Df_gse66863[,c(1,2,6)]
colnames(exp_g)[c(1,2)]<-c("gene","logFC_g")
colnames(methy)[2]<-"logFC_m"
exp_g<-exp_g%>%
  inner_join(methy,by="gene")
load("/Volumes/T7/Methylation/Download/04volcano.rda")
sub<-exp_g%>%inner_join(exp_methy_gene,by="Probe")
pdf(file='./04gse66863_volcano.pdf',width=4,height=4)
ggplot(exp_g, aes(x=logFC_m, y=logFC_g)) +
  geom_jitter(size=1,shape=21,fill="#c5c5c5") +
  theme_classic()+
  labs(x="DNA Methylation (log2FC)",y="Gene EXpression (log2FC)") +
  #scale_y_continuous(limits = c(0,yMax+0.2),expand = c(0,0)) +
  #scale_x_continuous(limits = c(-1.5,xMax)) +
  geom_hline(yintercept = 0, lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = 0, lty=4,col="grey",lwd=0.6) +
  geom_jitter(data=sub,aes(x=logFC_m, y=logFC_g, color=status),size=0.2,shape=21)+
  scale_color_manual(values=color)+
  #geom_text(data=sub,aes(x=logFC+0.1,y=-log10(adj.P.Val)+0.1,label=symbol),size=3)+
  theme(legend.position = "none",
        axis.text=element_text(size=12),
        axis.line=element_line(color="black"),
        axis.ticks=element_line(color="black"),
        axis.title=element_text(size=12,face="plain"),
        panel.background =element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank())
dev.off()
