rm(list=ls())
### Correlation Table S5-----
library(tidyverse)
load("/Users/huzixin/Data_Project//Methylation/00Probe_match.rda")
load("/Users/huzixin/Data_Project/Methylation/00Gene_exp_methy_match.rda")
TCGA_Methy<-data.table::fread("./TCGA_MyDMP.txt")
cis_probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe))
trans_probe<-unique(c(probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))
load("./02.TPM&Targets.Rda")
load("./04CHAMPline_import.Rda")
tcga_beta<-beta_norm
tcga_exp<-TPM_sym[,colnames(beta_norm)]
identical(colnames(tcga_exp),colnames(tcga_beta))
rm(beta_norm,Myload,TPM,TPM_sym)

probe<-cis_probe
probe_fl<-c()
esm<-c()
p.val<-c()
for(i in 1:length(probe)){
  a<-tcga_beta[probe[i],]
  gene<-TCGA_Methy[TCGA_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(tcga_exp)){
    b<-tcga_exp[gene,]%>%as.numeric()
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
probe_fl_tcga_cis<-probe_fl
esm_tcga_cis<-esm
pval_tcga_cis<-p.val
probe_fl_df_tcga_cis<-data.frame(Probe=probe_fl_tcga_cis,Estimate=esm_tcga_cis,P.val=pval_tcga_cis)
probe_fl_df_tcga_cis$cr<-rep("Cis",nrow(probe_fl_df_tcga_cis))

probe<-trans_probe
probe_fl<-c()
esm<-c()
p.val<-c()

for(i in 1:length(probe)){
  a<-tcga_beta[probe[i],]
  gene<-TCGA_Methy[TCGA_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(tcga_exp)){
    b<-tcga_exp[gene,]%>%as.numeric()
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

probe_fl_tcga_trans<-probe_fl
esm_tcga_trans<-esm
pval_tcga_trans<-p.val

probe_fl_df_tcga_trans<-data.frame(Probe=probe_fl_tcga_trans,
                                   Estimate=esm_tcga_trans,
                                   P.val=pval_tcga_trans)

probe_fl_df_tcga_trans$cr<-rep("Trans",nrow(probe_fl_df_tcga_trans))
probe_fl_tcga<-rbind(probe_fl_df_tcga_cis,probe_fl_df_tcga_trans)
colnames(probe_fl_tcga)[2:3]<-paste0(colnames(probe_fl_tcga)[2:3],"_tcga")

save(probe_fl_tcga,file="04Co_exp_methy_probe_fl_sp.rda")


# combine logFC-----
rm(list=ls())
load("./05Diffgene_tcga.rda")
table(Df_tcga$logFC>0)
table(Df_tcga$logFC<0)
TCGA_Methy<-data.table::fread("./TCGA_MyDMP.txt")
table(TCGA_Methy$logFC>0)
table(TCGA_Methy$logFC<0)
colnames(TCGA_Methy)
colnames(Df_tcga)
Methy<-TCGA_Methy%>%select("Probe","logFC","gene","feature","cgi","Enhancer","UCSC_CpG_Islands_Name")
Expr<-Df_tcga%>%select("Symbol","logFC","Chromosome","Start","End")
colnames(Expr)[c(2,1)]<-c("logFC_g","gene")
colnames(Methy)[2]<-"logFC_m"

tcga_Expr_Methy<-Expr%>%
  inner_join(Methy,by="gene")%>%
  dplyr::select("Probe","gene","logFC_g","logFC_m",
                "Chromosome","Start","End",
                "feature","cgi","Enhancer","UCSC_CpG_Islands_Name")
cis_upgene<-tcga_Expr_Methy[tcga_Expr_Methy$logFC_m>0&tcga_Expr_Methy$logFC_g>0,]
length(unique(cis_upgene$gene))
cis_downgene<-tcga_Expr_Methy[tcga_Expr_Methy$logFC_m<0&tcga_Expr_Methy$logFC_g<0,]
length(unique(cis_downgene$gene))
trans_downgene<-tcga_Expr_Methy[tcga_Expr_Methy$logFC_m>0&tcga_Expr_Methy$logFC_g<0,]
length(unique(trans_downgene$gene))
trans_upgene<-tcga_Expr_Methy[tcga_Expr_Methy$logFC_m<0&tcga_Expr_Methy$logFC_g>0,]
length(unique(trans_upgene$gene))
save(tcga_Expr_Methy,file="05exp_methy_logFC.rda")
#Gene----
library(tidyverse)
rm(list=ls())
load("./05Diffgene_tcga.rda")
TCGA_Methy<-data.table::fread("./TCGA_MyDMP.txt")
load("/Volumes/T7/Methylation/00Gene_exp_methy_match.rda")
colnames(TCGA_Methy)
colnames(Df_tcga) 
methy<-TCGA_Methy[,c(1,2,6,15)]
exp_g<-Df_tcga[,c(5,6,8,12)]
colnames(exp_g)[c(1,3)]<-c("gene","logFC_g")
colnames(methy)[2]<-"logFC_m"
exp_g<-exp_g%>%
  inner_join(methy,by="gene")


exp_methy_gene<-rbind(data.frame(gene=probe_exp_methy_d_cis$gene,
                                 Probe=probe_exp_methy_d_cis$Probe,
                                 status=rep("Cis_Hypo",length(probe_exp_methy_d_cis$gene))),
                      data.frame(gene=probe_exp_methy_d_trans$gene,
                                 Probe=probe_exp_methy_d_trans$Probe,
                                 status=rep("Trans_Hypo",length(probe_exp_methy_d_trans$gene))),
                      data.frame(gene=probe_exp_methy_u_cis$gene,
                                 Probe=probe_exp_methy_u_cis$Probe,
                                 status=rep("Cis_Hyper",length(probe_exp_methy_u_cis$gene))),
                      data.frame(gene=probe_exp_methy_u_trans$gene,
                                 Probe=probe_exp_methy_u_trans$Probe,
                                 status=rep("Trans_Hyper",length(probe_exp_methy_u_trans$gene))))
sub<-exp_g%>%inner_join(exp_methy_gene,by="Probe")
color<-c('#F8766D','#7cae00','#01bfc4','#c77cff')
save(exp_methy_gene,color,file="./04volcano.rda")
## Figure 3 F----
pdf(file='./04tcga_volcano.pdf',width=4,height=4)
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
scale_fill_manual(values=c('#F8766D','#7cae00','#01bfc4','#c77cff'),
                  labels=c("Trans-regulation,Hypo-methylation","Cis-regulation,Hyper-methylation",
                           "Cis-regulation,Hypo-methylation","Trans-regulation,Hyper-methylation"))