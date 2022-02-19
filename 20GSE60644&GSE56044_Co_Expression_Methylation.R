library(tidyverse)
library(ggpubr)
rm(list=ls())
load("/Users/huzixin/Data_Project//Methylation/00Gene_exp_methy_match.rda")
load("01Importdata.rda")
load("02Methylation.Rda")
GSE56044_Methy<-data.table::fread("./GSE56044_MyDMP.txt")
cis_probe<-unique(c(probe_exp_methy_d_cis$Probe,probe_exp_methy_u_cis$Probe))
trans_probe<-unique(c(probe_exp_methy_d_trans$Probe,probe_exp_methy_u_trans$Probe))

gse56044_beta<-myNorm
gse56044_exp<-norm[,colnames(myNorm)]
identical(colnames(gse56044_beta),colnames(gse56044_exp))
save(gse56044_beta,gse56044_exp,file="./04exp56044_methy60644_data")

rm(myNorm,norm,non_norm,pda_co,xxx)
ggdensity(gse56044_beta[probe[1],])
shapiro.test(gse56044_beta[probe[1],])
probe<-cis_probe
probe_fl<-c()
esm<-c()
p.val<-c()

for(i in 1:length(probe)){
  a<-gse56044_beta[probe[i],]
  gene<-GSE56044_Methy[GSE56044_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(gse56044_exp)){
    b<-gse56044_exp[gene,]%>%as.numeric()
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
probe_fl_56044_cis<-probe_fl
esm_56044_cis<-esm
pval_56044_cis<-p.val

probe_fl_df_56044_cis<-data.frame(Probe=probe_fl_56044_cis,
                                  Estimate=esm_56044_cis,
                                  P.val=pval_56044_cis)
probe_fl_df_56044_cis$cr<-rep("Cis",nrow(probe_fl_df_56044_cis))
probe<-trans_probe
probe_fl<-c()
esm<-c()
p.val<-c()

for(i in 1:length(probe)){
  a<-gse56044_beta[probe[i],]
  gene<-GSE56044_Methy[GSE56044_Methy$Probe==probe[i],"gene"]%>%as.character()
  if(str_detect(gene,"") & gene%in%rownames(gse56044_exp)){
    b<-gse56044_exp[gene,]%>%as.numeric()
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
probe_fl_56044_trans<-probe_fl
esm_56044_trans<-esm
pval_56044_trans<-p.val

probe_fl_df_56044_trans<-data.frame(Probe=probe_fl_56044_trans,
                                    Estimate=esm_56044_trans,
                                    P.val=pval_56044_trans)
probe_fl_df_56044_trans$cr<-rep("Trans",nrow(probe_fl_df_56044_trans))
probe_fl_56044<-rbind(probe_fl_df_56044_cis,probe_fl_df_56044_trans)
colnames(probe_fl_56044)[2:3]<-paste0(colnames(probe_fl_56044)[2:3],"_gse56044")
save(probe_fl_56044,file="04Co_exp_methy_probe_fl_sp.rda")

# combine logFC-----
rm(list=ls())
load("./05Diffgene_gse60644.rda")
table(Df_gse60644$up_down)
GSE56044_Methy<-data.table::fread("./GSE56044_MyDMP.txt")
colnames(GSE56044_Methy)
colnames(Df_gse60644)
Methy<-GSE56044_Methy%>%select("Probe","logFC","gene")
Expr<-Df_gse60644%>%select("Symbol","logFC")
colnames(Expr)[c(2,1)]<-c("logFC_g","gene")
colnames(Methy)[2]<-"logFC_m"

gse60644_Expr_Methy<-Expr%>%
  inner_join(Methy,by="gene")%>%
  dplyr::select("Probe","gene","logFC_g","logFC_m")

cis_genedown<-gse60644_Expr_Methy[gse60644_Expr_Methy$logFC_g<0&gse60644_Expr_Methy$logFC_m<0,]
cis_geneup<-gse60644_Expr_Methy[gse60644_Expr_Methy$logFC_g>0&gse60644_Expr_Methy$logFC_m>0,]
trans_genedown<-gse60644_Expr_Methy[gse60644_Expr_Methy$logFC_g<0&gse60644_Expr_Methy$logFC_m>0,]
trans_geneup<-gse60644_Expr_Methy[gse60644_Expr_Methy$logFC_g>0&gse60644_Expr_Methy$logFC_m<0,]
save(gse60644_Expr_Methy,file="05exp_methy_logFC.rda")

# Figure 3 H Volcano----
library(tidyverse)
rm(list=ls())
load("./05Diffgene_gse60644.rda")
Methy<-data.table::fread("./GSE56044_MyDMP.txt")

colnames(Methy)
colnames(Df_gse66863) 
methy<-Methy[,c(1,2,6,15)]
exp_g<-Df_gse60644[,c(1,2,6)]
colnames(exp_g)[c(1,2)]<-c("gene","logFC_g")
colnames(methy)[2]<-"logFC_m"
exp_g<-exp_g%>%
  inner_join(methy,by="gene")
load("/Volumes/T7/Methylation/Download/04volcano.rda")
sub<-exp_g%>%inner_join(exp_methy_gene,by="Probe")
pdf(file='./04gse56044_volcano.pdf',width=4,height=4)
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
