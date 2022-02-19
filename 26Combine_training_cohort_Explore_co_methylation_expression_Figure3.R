## Explore_coProbe&Gene----
#Figure 3 E I Number_Cis or Trans??----
load("./00Gene_exp_methy_match.rda")
exp_meth<-data.frame(Up=c(nrow(probe_exp_methy_d_trans),nrow(probe_exp_methy_u_cis)),
                     Down=c(nrow(probe_exp_methy_d_cis),nrow(probe_exp_methy_u_trans)))
rownames(exp_meth)<-c("Down","Up")
exp_meth<-as.table(as.matrix(exp_meth))
chisq.test(exp_meth)$expected
# G_up    G_down
# M_down 1446.3899 2359.6101
# M_up    395.6101  645.3899
# G_up G_down
# M_down 1519   2287
# M_up    323    718
#X-squared = 27, df = 1, p-value = 2.035e-07
exp_meth<-as.data.frame(exp_meth)%>%reshape::melt()
exp_meth<-exp_meth[,-3]
colnames(exp_meth)<-c("Methylation","Expression","Count")
library(ggpubr)
pdf("./Methy_Diff/01meth_exp_proportion.pdf",width=3,height=4)
ggbarplot(exp_meth, x = "Expression", y = "Count",
          fill = "Methylation", color = "Methylation", 
          palette = c("gray", "black"),
          label = TRUE, lab.col = "white", lab.pos = "in")
dev.off()
exp_meth$Status<-c("Trans-regulation,Hypo-methylation","Cis-regulation,Hyper-methylation",
                   "Cis-regulation,Hypo-methylation","Trans-regulation,Hyper-methylation")
exp_meth$share<-exp_meth$Count*100/sum(exp_meth$Count)
pdf("./Methy_Diff/01meth_exp_proportion_circle.pdf",width=12,height=4)
ggplot(exp_meth,aes('', y = Count, fill = Status))+
  geom_bar(width = 1, size = 1,stat = 'identity')+
  coord_polar(theta = 'y')+
  geom_text(aes(label = paste0(round(share,1),'%')),
            position = position_stack(vjust = .5),size=5)+
  labs(x = NULL, y = NULL, fill = NULL)+
  scale_fill_manual(values=c('#F8766D','#7cae00','#01bfc4','#c77cff'))+
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=12))
dev.off()
#Figure 3 L Chromosome---- 
load("./00Gene_exp_methy_match.rda")
colnames(probe_exp_methy_d_cis)
chr<-paste0("chr",as.character(1:22))
Feat_d_cis<-probe_exp_methy_d_cis[,c(1,2,9,12,13)]%>%
  add_count(Chromosome,name = "CHR_count")%>%
  dplyr::select("Chromosome","CHR_count")%>%
  distinct(Chromosome,.keep_all = T)
Feat_d_cis$status<-rep("Cis-regulation,Hypo-methylation",nrow(Feat_d_cis))

Feat_d_trans<-probe_exp_methy_d_trans[,c(1,2,9,12,13)]%>%
  add_count(Chromosome,name = "CHR_count")%>%
  dplyr::select("Chromosome","CHR_count")%>%
  distinct(Chromosome,.keep_all = T)
Feat_d_trans$status<-rep("Trans-regulation,Hypo-methylation",nrow(Feat_d_trans))

Feat_u_trans<-probe_exp_methy_u_trans[,c(1,2,9,12,13)]%>%
  add_count(Chromosome,name = "CHR_count")%>%
  dplyr::select("Chromosome","CHR_count")%>%
  distinct(Chromosome,.keep_all = T)
Feat_u_trans$status<-rep("Trans-regulation,Hyper-methylation",nrow(Feat_u_trans))

Feat_u_cis<-probe_exp_methy_u_cis[,c(1,2,9,12,13)]%>%
  add_count(Chromosome,name = "CHR_count")%>%
  dplyr::select("Chromosome","CHR_count")%>%
  distinct(Chromosome,.keep_all = T)
Feat_u_cis$status<-rep("Cis-regulation,Hyper-methylation",nrow(Feat_u_cis))

Feat<-rbind(Feat_d_cis,Feat_u_cis,Feat_d_trans,Feat_u_trans)
colnames(Feat)
pdf("./Methy_Diff/01chromosome.pdf",width=5,height=8)
ggdotchart(Feat, x = "Chromosome", y = "CHR_count",
           color = "status",                     
           add = "segments",                                 
           group = "Chromosome",                               
           dot.size = 6,   
           label = "CHR_count",  
           font.label = list(color = "white", size = 8, vjust = 0.5), 
           ggtheme = theme_pubr())+
  scale_fill_manual(values=c('#F8766D','#7cae00','#01bfc4','#c77cff'))+
  theme(axis.text.y  =  element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=8),
        axis.line.y = element_blank(),
        axis.title  =  element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position="none",
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(0, "cm"))+
  facet_wrap("status",nrow = 4)

dev.off()

#Figure 3 J Feature:Body TSS1500---- 
colnames(probe_exp_methy_d_cis)
Feat_d_cis<-probe_exp_methy_d_cis[,c(1,2,9,12)]%>%
  mutate(feature=ifelse(feature=="TSS200","TSS1500","non_TSS1500"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_d_cis$status<-rep("Cis-regulation,\nHypo-methylation",nrow(Feat_d_cis))
Feat_u_trans<-probe_exp_methy_u_trans[,c(1,2,9,12)]%>%
  mutate(feature=ifelse(feature=="TSS200","TSS1500","non_TSS1500"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_u_trans$status<-rep("Trans-regulation,\nHyper-methylation",nrow(Feat_u_trans))
Feat_d_trans<-probe_exp_methy_d_trans[,c(1,2,9,12)]%>%
  mutate(feature=ifelse(feature=="TSS200","TSS1500","non_TSS1500"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_d_trans$status<-rep("Trans-regulation,\nHypo-methylation",nrow(Feat_d_trans))
Feat_u_cis<-probe_exp_methy_u_cis[,c(1,2,9,12)]%>%
  mutate(feature=ifelse(feature=="TSS200","TSS1500","non_TSS1500"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_u_cis$status<-rep("Cis-regulation,\nHyper-methylation",nrow(Feat_u_cis))


Feat<-rbind(Feat_d_cis,Feat_u_cis,Feat_d_trans,Feat_u_trans)
Feat_wider<-Feat%>%
  pivot_wider(names_from = "status",values_from = "Feat_count")%>%
  column_to_rownames("feature")
Feat_wider<-as.table(as.matrix(as.data.frame(Feat_wider)))
##cis vs trans
#TSS1500 vs body 
x<-data.frame(c(375+8,1352+275),c(340+114,674+297))
#TSS1500 vs non-TS1500
x<-data.frame(c(375+8,1912+315),c(340+114,1179+604))
#Body vs None body
x<-data.frame(c(1352+275,935+48),c(674+297，845+421))
colnames(x)<-c("Cis","Trans")
rownames(x)<-c("TSS1500","Body")
rownames(x)<-c("TSS1500","non_TSS1500")
rownames(x)<-c("Body","non_Body")
chisq.test(x)$expected
# X-squared = 26.245, df = 1, p-value = 3.007e-07 Trans TSS1500
# X-squared = 73.494, df = 1, p-value < 2.2e-16
# X-squared = 21.283, df = 1, p-value = 3.963e-06 CIs Body

x<-data.frame(c(375+340,1352+674),c(8+114,275+297))
x<-data.frame(c(375+340,1912+1179),c(8+114,315+604))
x<-data.frame(c(1352+674,935+845),c(275+297,421+48))
colnames(x)<-c("Hypo","Hyper")
rownames(x)<-c("TSS1500","Body")
rownames(x)<-c("TSS1500","non_TSS1500")
rownames(x)<-c("Body","non_Body")
#X-squared = 258.05, df = 1, p-value < 2.2e-16
#X-squared = 28.08, df = 1, p-value = 1.164e-07
#X-squared = 0.89944, df = 1, p-value = 0.3429

Trans TSS1500
Cis Body
pdf("./Methy_Diff/01feature.pdf",width=6,height=4)
ggbarplot(Feat, x = "status", y = "Feat_count",
          fill = "feature", color="black",
          label = F, lab.col = "white", lab.pos = "in")+
  scale_fill_brewer(palette ="Paired")+
  guides(fill=guide_legend(title=""))+
  labs(y="Count")+
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(size=10),
    axis.ticks.x = element_blank(),
    axis.title.x =  element_blank())

dev.off()



#Figure 3 K CGI shore----
table(Feat_d_cis$Feat_cgi)
Feat_d_cis<-probe_exp_methy_d_cis[,c(1,2,12,13)]%>%
  filter(cgi%in%c("shore","island"))%>%
  mutate(cgi=ifelse(cgi=="shore","shore","non_shore"))%>%
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi_count","cgi")

Feat_d_cis$status<-rep("Cis-regulation,\nHypo-methylation",nrow(Feat_d_cis))
Feat_u_trans<-probe_exp_methy_u_trans[,c(1,2,12,13)]%>%
  filter(cgi%in%c("shore","island"))%>%
  mutate(cgi=ifelse(cgi=="shore","shore","non_shore"))%>%
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi_count","cgi")

Feat_u_trans$status<-rep("Trans-regulation,\nHyper-methylation",nrow(Feat_u_trans))
Feat_d_trans<-probe_exp_methy_d_trans[,c(1,2,12,13)]%>%
  filter(cgi%in%c("shore","island"))%>%
  mutate(cgi=ifelse(cgi=="shore","shore","non_shore"))%>%
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi_count","cgi")

Feat_d_trans$status<-rep("Trans-regulation,\nHypo-methylation",nrow(Feat_d_trans))

Feat_u_cis<-probe_exp_methy_u_cis[,c(1,2,12,13)]%>%
  filter(cgi%in%c("shore","island"))%>%
  mutate(cgi=ifelse(cgi=="shore","shore","non_shore"))%>%
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi_count","cgi")

Feat_u_cis$status<-rep("Cis-regulation,\nHyper-methylation",nrow(Feat_u_cis))

Feat<-rbind(Feat_d_cis,Feat_u_cis,Feat_d_trans,Feat_u_trans)

Feat_wider<-Feat%>%
  pivot_wider(names_from = "status",values_from = "cgi_count")%>%
  column_to_rownames("cgi")

#shore vs island hyper vs hypo
tab<-data.frame(shore=c(679+633,82+190),island=c(498+350,87+44))
rownames(tab)<-c("Hypo","Hyper")
#X-squared = 6.2781, df = 1, p-value =0.01222
# island had more hypo

tab<-data.frame(shore=c(82+633,679+190),island=c(498+87,350+44))
rownames(tab)<-c("cis","trans")
#shore had more trans
#X-squared = 51.129, df = 1, p-value =8.648e-13



tab<-data.frame(island=c(243,202),non=c(340-243,375-202))
tab<-data.frame(island=c(247,303),non=c(674-247,1352-303))
tab<-data.frame(island=c(679,633),non=c(1519-679,2287-633))

rownames(tab)<-c("Trans","Cis")
tab<-as.table(as.matrix(tab))
chisq.test(tab)$expected
chisq.test(tab)
# p-value = 0.008547
# island      non
# Trans 52.78322 287.2168
# Cis   58.21678 316.7832
# island non
# Trans     66 274
# Cis       45 330
#p-value = 1.826e-06
# shore non
# Trans    243  97
# Cis      202 173
# shore      non
# Trans 211.6084 128.3916
# Cis   233.3916 141.6084

pdf("./Methy_Diff/01meth_exp_tss1500cgi.pdf",width=6,height=4)
ggbarplot(Feat, x = "status", y = "cgi_count",
          fill = "cgi", color="black",
          label = F, lab.col = "white", lab.pos = "in")+
  scale_fill_brewer(palette ="Set3")+
  guides(fill=guide_legend(title=""))+
  labs(y="Count")+
  theme(
    panel.border = element_blank(),
    axis.text.x = element_text(size=10),
    axis.ticks.x = element_blank(),
    axis.title.x =  element_blank())
dev.off()

colnames(probe_exp_methy_d_cis)
Feat_d_cis<-probe_exp_methy_d_cis[,c(1,2,12,13)]%>%
  mutate(cgi=ifelse(cgi=="island","island","non_island") ) %>% 
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi","cgi_count")
Feat_d_cis$status<-rep("Cis-regulation,\nHypo-methylation",nrow(Feat_d_cis))

Feat_u_trans<-probe_exp_methy_u_trans[,c(1,2,12,13)]%>%
  mutate(cgi=ifelse(cgi=="island","island","non_island") ) %>% 
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi","cgi_count")

Feat_u_trans$status<-rep("Trans-regulation,\nHyper-methylation",nrow(Feat_u_trans))

Feat_d_trans<-probe_exp_methy_d_trans[,c(1,2,12,13)]%>%
  mutate(cgi=ifelse(cgi=="island","island","non_island") ) %>% 
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi","cgi_count")

Feat_d_trans$status<-rep("Trans-regulation,\nHypo-methylation",nrow(Feat_d_trans))

Feat_u_cis<-probe_exp_methy_u_cis[,c(1,2,12,13)]%>%
  mutate(cgi=ifelse(cgi=="island","island","non_island") ) %>% 
  add_count(cgi,name = "cgi_count")%>%
  distinct(cgi,.keep_all = T)%>%
  dplyr::select("cgi","cgi_count")

Feat_u_cis$status<-rep("Cis-regulation,\nHyper-methylation",nrow(Feat_u_cis))
Feat<-rbind(Feat_d_cis,Feat_u_cis,Feat_d_trans,Feat_u_trans)
sum(Feat$cgi_count[Feat$status%in%c("Trans-regulation,\nHypo-methylation","Cis-regulation,\nHypo-methylation")])#3806
Feat$cgi_count[(Feat$status%in%c("Trans-regulation,\nHypo-methylation","Cis-regulation,\nHypo-methylation"))& Feat$cgi=="island"]#350+498=848 
sum(Feat$cgi_count[Feat$status%in%c("Cis-regulation,\nHyper-methylation","Trans-regulation,\nHyper-methylation")])#1041
Feat$cgi_count[Feat$status%in%c("Cis-regulation,\nHyper-methylation","Trans-regulation,\nHyper-methylation")& Feat$cgi=="island"]#87+44=131
tab<-data.frame(island=c(848,131),non=c(3806-848,1041-131))
rownames(tab)<-c("hypo","hyper")
tab<-as.table(as.matrix(tab))
chisq.test(tab)$expected
chisq.test(tab)
#X-squared = 47.083, df = 1, p-value = 6.805e-12










## Promoter----
table(probe_exp_methy_d_cis$feature)
Feat_d_cis<-probe_exp_methy_d_cis[,c(1,2,9,12)]%>%
  filter(!(feature=="3'UTR")) %>% 
  mutate(feature=ifelse(feature=="Body","Body","Promoter"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_d_cis$status<-rep("Cis-regulation,\nHypo-methylation",nrow(Feat_d_cis))
Feat_u_trans<-probe_exp_methy_u_trans[,c(1,2,9,12)]%>%
  filter(!(feature=="3'UTR")) %>% 
  mutate(feature=ifelse(feature=="Body","Body","Promoter"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_u_trans$status<-rep("Trans-regulation,\nHyper-methylation",nrow(Feat_u_trans))
Feat_d_trans<-probe_exp_methy_d_trans[,c(1,2,9,12)]%>%
  filter(!(feature=="3'UTR")) %>% 
  mutate(feature=ifelse(feature=="Body","Body","Promoter"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_d_trans$status<-rep("Trans-regulation,\nHypo-methylation",nrow(Feat_d_trans))
Feat_u_cis<-probe_exp_methy_u_cis[,c(1,2,9,12)]%>%
  filter(!(feature=="3'UTR")) %>% 
  mutate(feature=ifelse(feature=="Body","Body","Promoter"))%>%
  add_count(feature,name = "Feat_count")%>%
  distinct(feature,.keep_all = T)%>%
  dplyr::select("feature","Feat_count")
Feat_u_cis$status<-rep("Cis-regulation,\nHyper-methylation",nrow(Feat_u_cis))


Feat<-rbind(Feat_d_cis,Feat_u_cis,Feat_d_trans,Feat_u_trans)
Feat_wider<-Feat%>%
  pivot_wider(names_from = "status",values_from = "Feat_count")%>%
  column_to_rownames("feature")

x<-data.frame(c(774+28,1352+275),c(750+390,674+297))
x<-data.frame(c(774+750,1352+674),c(28+390,275+297))
colnames(x)<-c("Cis","Trans")
colnames(x)<-c("Hypo","Hyper")
rownames(x)<-c("Promoter","Body")

chisq.test(x)$expected
# cis vs trans X-squared = 202.34, df = 1, p-value < 2.2e-16
# Hyper vs Hypo X-squared = 0.13066, df = 1, p-value = 0.7178


