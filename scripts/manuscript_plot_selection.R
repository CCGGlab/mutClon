# Load results
load("results/data/manuscript_selection.RData")
library(ggrepel)

#############
# Volcano 
#############

# Data frame
dNdS_volc_df<- data.frame(
  gene_id= rep(rownames(dNdS$PM$skin),2),
  organ=rep(c("skin","oral"),each=nrow(dNdS$PM$skin)),
  dNdS=c(dNdS$PM$skin$dNdS,dNdS$PM$oral$dNdS),
  p_dNdS=c(dNdS$PM$skin$p_dNdS,dNdS$PM$oral$p_dNdS), 
  q_dNdS=c(dNdS$PM$skin$q_dNdS,dNdS$PM$oral$q_dNdS),
  dNndS=c(dNdS$PM$skin$dNndS,dNdS$PM$oral$dNndS),
  p_dNndS=c(dNdS$PM$skin$p_dNndS,dNdS$PM$oral$p_dNndS), 
  q_dNndS=c(dNdS$PM$skin$q_dNndS,dNdS$PM$oral$q_dNndS),
  dPP2=c(PP2$PM$skin$obs-PP2$PM$skin$exp, PP2$PM$oral$obs-PP2$PM$oral$exp),
  p_PP2=c(PP2$PM$skin$p, PP2$PM$oral$p),
  q_PP2=c(PP2$PM$skin$q, PP2$PM$oral$q)
)
dNdS_volc_df<- dNdS_volc_df[grep("HLA.*",dNdS_volc_df$gene_id, invert = T),]

# Gene types seperately
dNdS_volc_df_type<- dNdS_volc_df[dNdS_volc_df$gene_id%in%c("all","driver","immune","housekeeping"),] 
dNdS_volc_df<- dNdS_volc_df[!dNdS_volc_df$gene_id%in%c("all","driver","immune","housekeeping"),] 

# Plot dN/dS
dNdS_volc_df$isSign<- dNdS_volc_df$q_dNdS<0.1
dNdS_volc_df$isLabel<- dNdS_volc_df$gene_id%in%unique(dNdS_volc_df$gene_id[dNdS_volc_df$isSign])

p_volc_dNdS<- ggplot(dNdS_volc_df, aes(x = dNdS, y = -log10(p_dNdS), color= isSign, key=gene_id)) +
  geom_point(size=2) +
  xlab("dN/dS") +
  ylab("-log10(P)") +
  scale_color_manual(values = c("#d0d3d4", "#3498db")) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    strip.text = element_text(size=8),
    axis.line = element_line(colour = "black")
  ) +
  facet_grid(.~organ) +
  geom_text_repel(data = subset(dNdS_volc_df, isLabel),
                  aes(label = gene_id, fontface=3),
                  size = 2,
                  colour="black",
                  box.padding   = 0.5,
                  point.padding = 0.5,
                  max.overlaps = 25,
                  segment.color = 'grey50')
p_volc_dNdS<- p_volc_dNdS +
  geom_point(data = subset(dNdS_volc_df_type, gene_id=="all"), aes(x = dNdS, y = -log10(p_dNdS), key=gene_id), colour="blue", shape = 3, cex=3)
p_volc_dNdS

# Plot PP2
dNdS_volc_df$isSign<- dNdS_volc_df$q_PP2<0.1 
dNdS_volc_df$isLabel<- dNdS_volc_df$gene_id%in%unique(dNdS_volc_df$gene_id[dNdS_volc_df$isSign])
p_volc_PP2<- ggplot(dNdS_volc_df, aes(x = dPP2, y = -log10(p_PP2), color= isSign, key=gene_id)) +
  geom_point() +
  scale_color_manual(values = c("#d0d3d4", "#3498db")) +
  xlab("delta PP2") +
  ylab("-log10(P)") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    strip.text = element_text(size=8),
    axis.line = element_line(colour = "black")
  ) +
  facet_grid(.~organ) +
  geom_text_repel(data = subset(dNdS_volc_df, isLabel),
                  aes(label = gene_id, fontface=3),
                  size = 2,
                  colour="black",
                  box.padding   = 0.5,
                  point.padding = 0.5,
                  max.overlaps = 25,
                  segment.color = 'grey50')
p_volc_PP2<- p_volc_PP2 +
  geom_point(data = subset(dNdS_volc_df_type, gene_id=="all"), aes(x = dPP2, y = -log10(p_PP2), key=gene_id), colour="blue", shape = 3, cex=3)
p_volc_PP2

#################
# dNdS barplots
#################

# Dataframe with results
genes<- na.omit(c(unique(dNdS_volc_df$gene_id[dNdS_volc_df$isSign]),"all","immune", "driver", "housekeeping"))
dNdS_df<- data.frame(
  gene = rep(factor(rep(genes,length(dNdS)),levels=genes[order(dNdS$PM$skin[genes,"p_dNdS"])]),3), # Order based on PM signfiicance in plot
  dNdS = c(as.numeric(sapply(dNdS, function(x) x$skin[genes,"dNdS"])),as.numeric(sapply(dNdS, function(x) x$skin[genes,"dNndS"])), as.numeric(sapply(dNdS, function(x) x$skin[genes,"dNmdS"]))),
  p_dNdS = c(as.numeric(sapply(dNdS, function(x) x$skin[genes,"p_dNdS"])), as.numeric(sapply(dNdS, function(x) x$skin[genes,"p_dNndS"])), as.numeric(sapply(dNdS, function(x) x$skin[genes,"p_dNmdS"]))),
  q_dNdS = c(as.numeric(sapply(dNdS, function(x) x$skin[genes,"q_dNdS"])), as.numeric(sapply(dNdS, function(x) x$skin[genes,"q_dNndS"])), as.numeric(sapply(dNdS, function(x) x$skin[genes,"q_dNmdS"]))),
  source = rep(factor(rep(c("PM","IM","SCC"),each=length(genes)),levels=c("PM","IM","SCC")),3),
  type = factor(rep(c("dN/dS","dNons/dS", "dMiss/dS"), each=length(dNdS)*length(genes)), levels=c("dN/dS", "dMiss/dS", "dNons/dS"))
)
dNdS_df$sign<- "NS"
dNdS_df$sign[dNdS_df$p_dNdS<0.05]<- "*"
dNdS_df$sign[dNdS_df$p_dNdS<0.01]<- "**"
dNdS_df$sign[dNdS_df$p_dNdS<0.001]<- "***"

dNdS_df$gene_type<- "individual"
for(type in c("all","immune", "driver", "housekeeping")) dNdS_df$gene_type[dNdS_df$gene==type]<- type
dNdS_df$gene_type<- factor(dNdS_df$gene_type, levels=c("all","individual","driver","immune","housekeeping"))
  
# Plot dN/dS by group, housekeeping not useful, only pseudocounts
p_dNdS_group <- ggplot(subset(dNdS_df,type!="dMiss/dS"&gene_type!="individual"&gene_type!="all"&source=="PM"), aes(x=gene, y=dNdS, fill=gene_type)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  xlab("") +
  ylab("dN/dS") +
  facet_grid(.~type) +
  scale_fill_manual(name="", values = c("#c06c84ff","#355c7dff","#f8b195ff")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7),
    strip.placement = "outside",
    strip.text = element_text(size=7),
    strip.background = element_blank(),
    axis.line.y = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_text(
    aes(x = gene, y = dNdS, group=source, label = format(p_dNdS, digits = 2, scientific = T)),
    position = position_dodge(0.9),
    vjust=-0.5, size=2
  ) +
  geom_hline(yintercept = 1, linetype="dashed")
p_dNdS_group 

# Plot dN/dS (sorted following p values PM study)
p_dNdS <- ggplot(subset(dNdS_df,type!="dMiss/dS"&gene_type=="individual"), aes(x=gene, y=dNdS, fill=source)) +
  geom_bar(stat="identity", color="white", position=position_dodge())+
  xlab("") +
  ylab("dN/dS") +
  facet_grid(type~., scales="free", switch = "y") +
  # scale_fill_manual(name="", values = c("blue","grey","red")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size=8, face = "italic"),
    axis.text.y = element_text(size=6),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size=7),
    strip.background.y = element_blank(),
    axis.line.y = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_text(
      aes(x = gene, y = dNdS, group=source, label = sign),
      position = position_dodge(0.9),
      vjust=-0.5, size=2
    ) +
  geom_hline(yintercept = 1, linetype="dashed")
p_dNdS 


#################
# PP2 violin plots
#################

# Simulations
load("data/PM_simulations.rds")

# Compare observed
PP2_obs_df<- as.data.frame(maf[,c("gene","PP2","study")])
PP2_df<- PP2_obs_df
for(gene in genes){
  PP2_exp_df<- data.frame(
    gene=gene,
    PP2=sample(PP2_exp[[gene]]$SKCM,1000),
    study="Exp"
  )
  PP2_df<- rbind(PP2_df, PP2_exp_df)
}
PP2_df$study<- factor(PP2_df$study, levels=c("Exp", "PM", "IM2015", "SCC"))

library(ggpubr)
PP2_df_tmp<- PP2_df
PP2_df_tmp$PP2[PP2_df_tmp$study=="Exp"]<- NA
PP2_df_tmp_Exp<- PP2_df
PP2_df_tmp_Exp$PP2[PP2_df_tmp_Exp$study!="Exp"]<- NA

# genes<- genes[grep("immune|driver|housekeeping",genes, invert = T)]
# PP2_df_tmp$gene<- factor(PP2_df_tmp$gene, levels=c("NOTCH1", "TP53", "FAT1", "NOTCH2", "driver", "immune", "housekeeping"))
PP2_df_tmp$gene<- factor(PP2_df_tmp$gene, levels=levels(dNdS_df$gene))
genes<- genes[grep("immune|driver|housekeeping",genes, invert = T)]
p_PP2<- ggplot(subset(PP2_df_tmp, gene%in%genes), aes(x=study, y=PP2, color=study)) + 
  geom_violin(data=subset(PP2_df_tmp_Exp, gene%in%genes), mapping=aes(x=study, y=PP2), fill="grey", color="black") +
  geom_jitter(draw_quantiles = T, size=1, alpha=0.5) +
  stat_summary(data=subset(PP2_df, gene%in%genes), mapping=aes(x=study, y=PP2, color=study), fun.y=median, geom="crossbar", color="black", size = 0.15) +
  stat_compare_means(data=subset(PP2_df, gene%in%genes), mapping=aes(x=study, y=PP2, color=study), label = "p.signif", method = "wilcox.test",
                     ref.group = "Exp", label.y = 1, size=2) +    
  facet_grid(.~gene, switch = "x") +
  xlab("") +
  ylab("PolyPhen-2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.09)) +
  scale_color_manual(values = c("white","#F8766D","#00BA38","#619CFF")) +
  # ylim(0,1) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=6),
    axis.title = element_text(size=7),
    axis.line.y = element_line(colour = "black"),
    strip.text = element_text(size=8, face="italic"),
    strip.background = element_blank()  
    )
p_PP2

genes<- c("driver", "immune", "housekeeping")
p_PP2_group<- ggplot(subset(PP2_df_tmp, gene%in%genes&study%in%c("Exp","PM")), aes(x=study, y=PP2)) + 
  geom_violin(data=subset(PP2_df_tmp_Exp, gene%in%genes&study==c("Exp","PM")), mapping=aes(x=study, y=PP2), fill="grey", color="black") +
  geom_jitter(draw_quantiles = T, size=1, alpha=0.5, col="#c06c84ff") +
  stat_summary(data=subset(PP2_df, gene%in%genes&study==c("Exp","PM")), mapping=aes(x=study, y=PP2, color=study), fun.y=median, geom="crossbar", color="black", size = 0.15) +
  stat_compare_means(data=subset(PP2_df, gene%in%genes&study==c("Exp","PM")), mapping=aes(x=study, y=PP2, color=study), label = "p.format", method = "wilcox.test",
                     label.x = 1.5, label.y = 1.05, size=2) +    
  facet_grid(.~gene, switch = "x") +
  xlab("") +
  ylab("PolyPhen-2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.09)) +
  # scale_color_manual(values = c("white","#F8766D","#00BA38","#619CFF")) +
  # scale_fill_manual(name="", values = c("#c06c84ff","#355c7dff","#f8b195ff")) +
  # ylim(0,1) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size=7),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=6),
    axis.title = element_text(size=7),
    axis.line.y = element_line(colour = "black"),
    strip.text = element_text(size=8, face="italic"),
    strip.background = element_blank()  
  )
p_PP2_group

#############################
# Some numbers
#############################

# Overall values?

# All genes in skin
format(dNdS$PM$skin["all","n"]/dNdS$PM$skin["all","s"], digits=3)
# n/s = 3.68

format(dNdS$PM$skin["all","NS"], digits=3)
# 2.12

format(dNdS$PM$skin[c("all","driver","immune","housekeeping"),c("dNdS", "p_dNdS")], digits=3)
#         dNdS   p_dNdS
# all          1.725 3.69e-08
# driver       1.775 2.37e-08
# immune       0.814 7.23e-01
# housekeeping 0.648 8.14e-01

format(dNdS$PM$skin[c("all","driver","immune","housekeeping"),c("dNndS", "p_dNndS")], digits=3)
# dNndS  p_dNndS
# all           3.04 3.87e-09
# driver        3.30 7.06e-10
# immune        1.63 6.80e-01
# housekeeping  1.75 1.00e+00

format(dNdS$PM$skin[c("all","driver","immune","housekeeping"),c("dNmdS", "p_dNmdS")], digits=3)
# dNmdS  p_dNmdS
# all          1.798 6.26e-09
# driver       1.664 1.13e-06
# immune       0.813 7.28e-01
# housekeeping 0.700 7.85e-01

format(dNdS$PM$oral["all",c("dNdS", "p_dNdS")], digits=3)
# dNdS p_dNdS
# all 2.13 0.000133

format(dNdS$PM$skin[c("NOTCH1", "TP53", "FAT1", "NOTCH2", "BCORL1","CDKN2A","AJUBA"),c("dNdS", "q_dNdS", "dNndS","q_dNndS")], digits=3)
# dNdS   q_dNdS dNndS  q_dNndS
# NOTCH1 11.25 1.97e-07 43.28 1.76e-08
# TP53    8.38 9.01e-04 29.95 2.79e-04
# FAT1   10.41 8.39e-03 64.04 3.82e-05
# NOTCH2  1.64 1.00e+00  4.51 1.00e+00
# BCORL1  2.01 1.00e+00 28.17 2.21e-01
# CDKN2A  1.60 1.00e+00 18.83 3.51e-01
# AJUBA   1.44 1.00e+00 50.95 8.61e-02

# PP2 
genes<- c("NOTCH1", "TP53", "FAT1", "NOTCH2", "BCORL1", "CDKN2A","AJUBA","all","driver", "immune", "housekeeping")
PP2$PM$skin[genes,]
#              obs   exp            p            q
# NOTCH1       0.9990 0.088 1.753871e-15 1.596022e-13
# TP53         1.0000 0.102 2.852544e-11 1.297908e-09
# FAT1         0.9985 0.329 2.344791e-05 4.267519e-04
# NOTCH2       0.9990 0.326 1.265680e-02 1.279743e-01
# BCORL1       1.0000 0.275 3.398883e-03 5.154972e-02
# CDKN2A       0.9985 0.199 4.073295e-03 5.295283e-02
# AJUBA        1.0000 0.010 7.217681e-03 8.210113e-02
# all          0.4975 0.107 2.250334e-09 5.119511e-08
# driver       0.6020 0.120 1.468735e-10 4.455163e-09
# immune       0.0095 0.035 8.269520e-01 9.622831e-01
# housekeeping 0.0120 0.056 8.591195e-01 9.622831e-01

PP2$IM2015$skin[genes,]
# obs   exp            p            q
# NOTCH1       0.9990 0.088 2.786778e-61 1.950744e-59
# TP53         0.9990 0.102 7.055759e-13 9.878062e-12
# FAT1         0.9870 0.329 2.079971e-07 2.426632e-06
# NOTCH2       0.9990 0.326 2.280327e-23 3.990573e-22
# BCORL1           NA    NA           NA           NA
# CDKN2A       0.4955 0.199 5.444414e-01 7.249940e-01
# AJUBA        0.1900 0.010 2.065576e-01 4.606403e-01
# driver       0.6200 0.120 3.078729e-46 7.183701e-45
# immune           NA    NA           NA           NA
# housekeeping     NA    NA           NA           NA

PP2$SCC$skin[genes,]
# obs   exp            p            q
# NOTCH1       0.9990 0.088 2.815614e-13 9.925039e-12
# TP53         1.0000 0.102 1.952236e-17 9.175511e-16
# FAT1         0.9980 0.329 8.879268e-04 1.251977e-02
# NOTCH2       0.9965 0.326 4.280801e-05 8.622757e-04
# BCORL1           NA    NA           NA           NA
# CDKN2A       1.0000 0.199 5.978854e-13 1.686037e-11
# AJUBA        0.9980 0.010 4.778449e-02 2.807339e-01
# driver       0.5160 0.120 3.320867e-19 4.682422e-17
# immune       0.0740 0.035 1.358897e-01 5.029226e-01
# housekeeping 0.1760 0.056 1.113108e-01 4.373049e-01

####################
# Merge & save plots
####################
p_volc<- plot_grid(
  p_volc_dNdS, NULL, p_volc_PP2, NULL, NULL, NULL,
  ncol = 2,
  rel_widths = c(1,1),
  rel_heights = c(1,1,3),
  axis = "rlbt",
  align = "hv"
)
p_volc
ggsave("results/figs/mutClon_selection_volcano.pdf", p_volc, width = 178, height = 265, units = "mm")

p_driver_immune<- plot_grid(
  p_dNdS_group + theme(axis.title.y = element_blank()), p_PP2_group  + theme(axis.title.y = element_blank()),  
  ncol=1,
  align="v"
)
p_driver_immune
ggsave("results/figs/mutClon_selection_driver_immune.pdf", p_driver_immune, width = 178/3, height = 2*265/5, units = "mm")

p_genes<- plot_grid(
  p_dNdS,  
  p_PP2,
  NA,
  ncol=1,
  rel_heights = c(3,2,5),
  align="vh",
  axis="l",
  labels = c("a","b"),
  label_size = 8
)
p_genes
ggsave("results/figs/mutClon_selection_genes.pdf", p_genes, width = 178, height = 265, units = "mm")
