############################################################################
# Generate clonal map + compare size, frequency/cm2 and % surface with IM 
############################################################################

# Load
#######
mutClon_maf<- readRDS("data/mutClon_maf.rds")
source("scripts/functions/getMutClonMapCoordinates.R")
source("scripts/functions/createMutClonMap.R")

# Filter
##########

# Only focus on samples with biopsy sizes of 5mm (punch biopsies)
mutClon_maf<- mutClon_maf[mutClon_maf$biopsy=="P",]

# Merge mutations that occur on same samples
clone_size_mm2<- tapply(mutClon_maf$clone_size_mm2, mutClon_maf$mut_id, "sum")
mutClon_maf$clone_size_mm2<- clone_size_mm2[mutClon_maf$mut_id]
mutClon_maf<- mutClon_maf[!duplicated(mutClon_maf$mut_id)|is.na(mutClon_maf$mut_id),]

# Only focus on skin samples
mutClon_maf<- mutClon_maf[mutClon_maf$organ=="skin",]

# Clonal map
############

# Get cummulative biopsy size
sample_surf_t<- table(mutClon_maf$sample_alias, mutClon_maf$surface_cm2)>0
total_biopsy_size_PM<- 100*sum(as.numeric(colnames(sample_surf_t)) * t(sample_surf_t)) # in mm2

# Plot
maf_co<- getMutClonMapCoordinates(maf = mutClon_maf, total_biopsy_size = total_biopsy_size_PM, avoid_partial_overlaps = T, show_plot = T)
maf_co$gene_name<- as.character(maf_co$gene)
saveRDS(maf_co, "temp/maf_co.rds")

# genes_to_plot<- c("NOTCH1", "TP53", "FAT1", "NOTCH2")
genes_to_plot<- c("NOTCH1", "TP53", "FAT1", "NOTCH2", "BCORL1", "CDKN2A", "AJUBA")
maf_co$gene_name[!maf_co$gene_name%in%genes_to_plot]<- "Other"

geneCols<- c("#C3A0D3FF","#D0A9ADFF","#F9D0A4FF","#ED99A4FF","#B0C1D2FF",rgb(1,0.5,0,0.5),rgb(0,1,0,0.5),NA)
names(geneCols)<- c(genes_to_plot,"Other")
geneBorderCols<- c("#5C2685FF","#882E36FF","#F38928FF","#CE2125FF","#3A618FFF",rgb(1,0.5,0),rgb(0,1,0),"grey")
names(geneBorderCols)<- c(genes_to_plot,"Other")

# # Don't override, maf_co = random!
# pdf("results/figs/mutClon_skin_clonal_map.pdf")
# createMutClonMap(maf = maf_co, geneCols = geneCols, geneBorderCols = geneBorderCols, meta_biopsy_r = 10) # size of 10mm
# dev.off()

# Compare
############

# Load IM to compare
maf_IM2015<- readRDS("data/IM2015.rds")
load("~/data/IM2015_biopsy_info.RData")

# Merge
clone_df<- data.frame(
  gene = c(as.character(mutClon_maf$gene), maf_IM2015$gene_name),
  clone_size = c(mutClon_maf$clone_size, maf_IM2015$clone_size_pred),
  source = c(rep("PM",nrow(mutClon_maf)), rep("IM2015",nrow(maf_IM2015)))
)
clone_df<- clone_df[clone_df$gene%in%genes_to_plot,]
clone_df$gene<- factor(clone_df$gene, levels=(genes_to_plot))
clone_df$source<- factor(clone_df$source, levels=c("PM","IM2015"))

# Plot clone size
# Note that hardly any synonymous muts, so hard to compare s to n!
p_clone_size <- ggplot(data=clone_df, aes(x=gene, y=clone_size, fill=source)) +
  geom_boxplot(lwd=0.20, outlier.size = 0.15, position = position_dodge(preserve = "single")) +
  xlab("") +
  ylab(bquote("Clone size (" ~mm^2~ ")")) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.text.x=element_text(face="italic"),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1.2)) # Axis at 0
p_clone_size

# Plot % clones
clone_size_total<- tapply(clone_df$clone_size, list(clone_df$gene, clone_df$source), "sum") 
clon_prop<- 100*(clone_size_total/cbind(rep(sum(total_biopsy_size_PM),nrow(clone_size_total)),rep(sum(total_biopsy_size),nrow(clone_size_total))))[genes_to_plot,]
clon_prop_df<- melt(clon_prop, varnames = c("gene","source"))

p_prop <- ggplot(data=clon_prop_df, aes(x=gene, y=value, fill=source)) +
  geom_bar(stat="identity", color=NA, position=position_dodge())+
  xlab("") +
  ylab("% skin") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.text.x=element_text(face="italic"),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = c(0, 0)) # Axis at 0
p_prop

# Plot clone frequency
clone_n<- table(clone_df$gene, clone_df$source) 
clone_freq<- 100*(clone_n/cbind(rep(sum(total_biopsy_size_PM),nrow(clone_n)),rep(sum(total_biopsy_size),nrow(clone_n))))[genes_to_plot,]
clone_freq_df<- melt(clone_freq, varnames = c("gene","source"))

p_freq <- ggplot(data=clone_freq_df, aes(x=gene, y=value, fill=source)) +
  geom_bar(stat="identity", color=NA, position=position_dodge())+
  xlab("") +
  ylab(bquote("Clones per " ~cm^2)) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.line.x = element_blank(),
    axis.text.x=element_text(face="italic"),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    axis.ticks.x = element_blank()
  ) +  
  scale_y_continuous(expand = c(0, 0)) # Axis at 0
p_freq

# Numbers
##########

# Clones per cm2?
round(colSums(clone_freq),0)
# PM IM2015 
# 78    193 

round(clone_freq,0)
#         PM IM2015
# NOTCH1 35     97
# TP53   16     20
# FAT1   15     32
# NOTCH2  6     42
# BCORL1  2      0
# CDKN2A  2      1
# AJUBA   2      2

# Sizes? (in mm2)
round(tapply(clone_df$clone_size, list(clone_df$gene, clone_df$source), "median"),2)
#       PM IM2015
# NOTCH1 0.30   0.14
# TP53   0.33   0.18
# FAT1   0.23   0.09
# NOTCH2 0.30   0.10
# BCORL1 0.49     NA
# CDKN2A 0.32   0.10
# AJUBA  0.37   0.13

round(tapply(clone_df$clone_size, list(clone_df$gene, clone_df$source), "max"),2)
# PM    IM2015
# NOTCH1 1.90   2.06
# TP53   0.66   2.27
# FAT1   1.19   1.06
# NOTCH2 0.56   0.99
# BCORL1 0.68     NA
# CDKN2A 0.59   0.11
# AJUBA  0.64   0.40

clon_prop_df
# gene source     value
# 1  NOTCH1     PM 15.20888889
# 2    TP53     PM  5.72200000
# 3    FAT1     PM  4.68355556
# 4  NOTCH2     PM  2.05800000
# 5  BCORL1     PM  0.89111111
# 6  CDKN2A     PM  0.83666667
# 7   AJUBA     PM  0.72422222
# 8  NOTCH1 IM2015 22.61727945
# 9    TP53 IM2015  5.07067444
# 10   FAT1 IM2015  5.00215380
# 11 NOTCH2 IM2015  6.77466605
# 12 BCORL1 IM2015          NA
# 13 CDKN2A IM2015  0.05271746
# 14  AJUBA IM2015  0.33342171

# Merge plots
##############
library(cowplot)
p<- plot_grid(
  NULL, p_clone_size + theme(legend.position = "none", axis.text.x = element_blank()), NULL, p_prop+ theme(legend.position = "none", axis.text.x = element_blank()), NULL, p_freq + theme(legend.position = "none"),
  ncol=2,
  rel_widths = c(2,1),
  axis = "rlbt",
  # labels = c(NA,"B",NA,"C",NA,"D"),
  # label_size = 12,
  align = "hv"
)
p

ggsave("results/figs/mutClon_clonality.pdf", p, width = 2*178/3, height = 265/4, units = "mm")
