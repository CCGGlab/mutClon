##################################
# Basic description of TS results
##################################

# Load
#######
library(ggplot2)
library(cowplot)
mutClon_maf<- readRDS("data/mutClon_maf.rds")

# N samples?
#############

sample_location_t<- table(mutClon_maf$sample, mutClon_maf$location)

# Total?
nrow(sample_location_t) # 17

# Per location?
colSums(sample_location_t>0)
# cheek  eyelid nose 
# 5       5     7 

# Mutation calls
################

# Total?
nrow(mutClon_maf) #910

# How many in multiple samples? 24
mut_id_multi<- mutClon_maf$mut_id[duplicated(mutClon_maf$mut_id)]
# View(mutClon_maf[mutClon_maf$mut_id%in%mut_id_multi,])

# Per sample?
sample_t<- table(mutClon_maf$sample)
summary(as.numeric(sample_t))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11.00   40.00   52.00   53.53   65.00   92.00 

# VAF?
# tapply(mutClon_maf$VAF, mutClon_maf$sample, "summary")
summary(mutClon_maf$VAF)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00232 0.00570 0.00779 0.01077 0.01250 0.10600 

p_vaf<- ggplot(mutClon_maf, aes(x=VAF)) +
  geom_histogram(binwidth=0.001, color="black", fill="grey") +
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(0,0.05) +
  xlab("Variant allele frequency") +
  ylab("Mutation counts") +
  theme(
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black")
  )
p_vaf

# Mut loads
sample_t<- table(mutClon_maf$sample_alias)
sample_t
# 1E  2N  3N  4N  5O  6O  7E  8E  9E 10E 11N 12N 13N 14N 15O 16O 17O 
# 46  61  52  40  11  13  56  45  92  38  65  49  91  64  19  77  91 

sample_t<- table(mutClon_maf$sample)
tapply(as.numeric(sample_t), substr(names(sample_t),7,8), "median")
# 01 02 04 
# 46 61 19  
tapply(as.numeric(sample_t), list(substr(names(sample_t),1,4), substr(names(sample_t),7,8)), "median")
#      01 02   04
# PM01 46 52  12
# PM02 50.5 64.5 77

#######################################################
# Oncoplot
#######################################################

# Dataframe with mutation numbers per sample 
mut_sample_table<- table(mutClon_maf$gene, mutClon_maf$sample_alias) # Solution gives bad facet grid visualization
mut_sample_table<- mut_sample_table[order(rowSums(mut_sample_table),decreasing = T),] # Sort, to get top 20
genes_top10<- rownames(mut_sample_table)[1:20]
mut_sample_table<- mut_sample_table[order(rowSums(mut_sample_table),decreasing = F),] # Sort, to get correct order in geom tile
mut_sample_table_df<- as.data.frame(mut_sample_table)
colnames(mut_sample_table_df)<- c("gene", "sample", "n_mut")
mut_sample_table_df$col<- mut_sample_table_df$n_mut
mut_sample_table_df$col[mut_sample_table_df$col>5]<- 5 
mut_sample_table_df$n_mut[mut_sample_table_df$n_mut==0]<- NA
mut_sample_table_df$text_col<- "black"
mut_sample_table_df$text_col[mut_sample_table_df$col>=4]<- "white"
mut_sample_table_df$patient<- NA
for(pt in unique(mutClon_maf$patient)) mut_sample_table_df$patient[mut_sample_table_df$sample%in%as.character(unique(mutClon_maf$sample_alias[mutClon_maf$patient==pt]))]<- pt
mut_sample_table_df<- mut_sample_table_df[!is.na(mut_sample_table_df$patient),]
mut_sample_table_df$gene_type<- NA
for(t in c("driver","immune","housekeeping")) mut_sample_table_df$gene_type[mut_sample_table_df$gene%in%as.character(unique(mutClon_maf$gene[mutClon_maf$gene_type==t]))]<- t

p_mut<- ggplot(subset(mut_sample_table_df, gene%in%genes_top10), aes(sample, gene, fill= col)) + 
  geom_tile(colour = "white") +
  geom_text(aes(label = n_mut, colour=text_col), size=7*0.35) +
  scale_color_manual(breaks = c("black", "white"),
                     values=c("black", "white"), guide=guide_none()) +
  scale_fill_gradient(low="white", high="blue") +
  scale_x_discrete(position = "top") +
  xlab("") +
  ylab("") +
  labs(fill = "# mutations") +
  facet_grid(.~ patient, space = 'free_x', scales = 'free_x') +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background.x = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    axis.text.x = element_text(size=8),  
    axis.text.y = element_text(size=8, face="italic"),  
  )
p_mut

p_mut_all<- ggplot(subset(mut_sample_table_df, !is.na(n_mut)), aes(sample, gene, fill= col)) + 
  geom_tile(colour = "white") +
  geom_text(aes(label = n_mut, colour=text_col), size=7*0.35) +
  scale_color_manual(breaks = c("black", "white"),
                     values=c("black", "white"), guide=guide_none()) +
  scale_fill_gradient(low="white", high="blue") +
  scale_x_discrete(position = "top") +
  xlab("") +
  ylab("") +
  labs(fill = "# mutations") +
  facet_grid(gene_type ~ patient, space = 'free', scales = 'free') +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background.x = element_blank(),
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    axis.text.x = element_text(size=8),  
    axis.text.y = element_text(size=6, face="italic")
  )
p_mut_all

# Frequency of top 10
top10<- sort(table(mutClon_maf$gene), decreasing = T)[1:10]
top10
# NOTCH1   TP53  MUC17   FAT1   APOB NOTCH2 SPHKAP   FAT4  PREX2 TRIOBP 
# 113     52     42     37     35     28     25     24     24     20 
top10_norm<- top10/length(levels(mutClon_maf$sample_alias))
round(top10_norm,1)
# NOTCH1   TP53  MUC17   FAT1   APOB NOTCH2 SPHKAP   FAT4  PREX2 TRIOBP 
# 6.6    3.1    2.5    2.2    2.1    1.6    1.5    1.4    1.4    1.2 

# NOTCH1
summary(mut_sample_table_df[mut_sample_table_df$gene=="NOTCH1","n_mut"])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 1.000   4.000   6.000   6.647   9.000  16.000 

# Hom many genes mutated?
gene_t<- table(mutClon_maf$gene)
gene_t<- gene_t[gene_t!=0]
length(gene_t) # 119

# Examples
genes_eg<- gene_t[c("FAT4", "ERBB4", "NOTCH3", "BCORL1", "EGFR", "ERBB2")]
genes_eg/length(levels(mutClon_maf$sample_alias))
# FAT4     ERBB4    NOTCH3    BCORL1      EGFR     ERBB2 
# 1.4117647 0.9411765 0.8235294 0.3529412 0.2941176 0.4117647 

# Immune genes?
immune_genes<- unique(as.character(mutClon_maf$gene[mutClon_maf$gene_type=="immune"]))
genes_eg_i<- gene_t[immune_genes]
# sort(genes_eg_i/length(levels(mutClon_maf$sample_alias)), decreasing = T)
sort(genes_eg_i, decreasing = T)
# PSD4    MYO1E     NOS2   IL17RA     LAG3  ALOX15B  ALOX12B     CANX   IL15RA     CD86 
# 5        4        4        3        3        3        3        3        2        2 
# PIP5K1A     TLR3 PDCD1LG2  ARL14EP    CD274  TNFRSF9     CD80 TRAF3IP2   PIK3R2    ERAP1 
# 2        2        1        1        1        1        1        1        1        1 
# TAP2     CALR 
# 1        1 

# Add mutation burden & subst types on top
###########################################

# SUbst type cols
substTypes<-c("C>A","C>G","C>T","T>A","T>C","T>G")
substTypes_cols<-c("blue","black","red","grey","green","pink")
names(substTypes_cols)<-substTypes

# Sample subst df
sample_subst_table<- table(mutClon_maf$sample_alias, mutClon_maf$subst_type)
sample_subst_table_df<- as.data.frame(sample_subst_table)
colnames(sample_subst_table_df)<- c("sample", "subst", "n_mut")
sample_subst_table_df$subst<- factor(sample_subst_table_df$subst, levels=rev(substTypes))
for(pt in c("PM01","PM02")) sample_subst_table_df$patient[sample_subst_table_df$sample%in%as.character(unique(mutClon_maf$sample_alias[mutClon_maf$patient==pt]))]<- pt

# Stacked
p_subst<- ggplot(sample_subst_table_df, aes(x=sample, y=n_mut, fill=subst)) + 
  geom_bar(position="fill", stat="identity", ) +
  scale_fill_manual(breaks = substTypes, values=substTypes_cols) +
  facet_grid(.~ patient, space = 'free_x', scales = 'free_x', switch = 'x') + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("") +
  ylab("") +
  theme(
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.placement = 'outside',
    strip.background.x = element_blank(),
    legend.text = element_text(size=7),
    axis.text.x = element_text(size=8),  
    axis.text.y = element_text(size=6),  
    axis.line.y = element_line(colour = "black")
  )
p_subst

# Add variants & mut frequency on the right
###########################################

# Merge non point mutations to "other" group
mutClon_maf$variant[is.na(mutClon_maf$variant)|mutClon_maf$variant%in%c("frameshift deletion", "frameshift substitution")]<- "frameshift indel" 
mutClon_maf$variant[mutClon_maf$variant%in%c("nonframeshift deletion", "nonframeshift substitution", "startloss", "stoploss")]<- "other" 

# Variant cols
vars<- c("synonymous SNV", "nonsynonymous SNV", "stopgain", "frameshift indel", "other","DNV")
var_cols<- c("grey","#00A5CF","#DE1A1A", "#FFBF00", "black", "#29BF12")
names(vars)<- var_cols

gene_variant_table<- table(mutClon_maf$gene, mutClon_maf$variant)

# Filter top 20
gene_variant_table<- gene_variant_table[rev(genes_top10),]

# Gene Variant df
gene_variant_table_df<- as.data.frame(gene_variant_table)
colnames(gene_variant_table_df)<- c("gene", "variant", "n_mut")
gene_variant_table_df$variant<- as.character(gene_variant_table_df$variant)
gene_variant_table_df$variant[gene_variant_table_df$variant%in%c("stoploss", "startloss")]<- "stoploss/startloss"
gene_variant_table_df$variant<- factor(gene_variant_table_df$variant, levels=rev(vars))

# plot
p_var<- ggplot(gene_variant_table_df, aes(x=n_mut, y=gene, fill=variant)) + 
  geom_bar(position="stack", stat="identity", ) +
  scale_fill_manual(breaks = vars, values=var_cols) +
  xlab("# mutations") +
  scale_x_continuous(position = "top", expand=c(0,0)) +
  ylab("") +
  theme(
    legend.position = c(1, 0),
    legend.justification = c("right", "bottom"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.text = element_text(size=7),
    axis.text.x = element_text(size=6),
    axis.title.x = element_text(size=7),
    axis.text.y = element_text(size=8, face="italic") ,
    axis.line.x = element_line(colour = "black")
  )
p_var

# Mut burden
#############
sample_pt_t<- table(mutClon_maf$sample_alias, mutClon_maf$patient, mutClon_maf$gene_type)
sample_pt_df<- as.data.frame(sample_pt_t)
colnames(sample_pt_df)<- c("sample", "patient", "type", "mut_burden")
sample_pt_df<- sample_pt_df[sample_pt_df$mut_burden!=0,,]  

p_burden<- ggplot(sample_pt_df, aes(x=sample,y=mut_burden, fill=type)) +
  # geom_bar(stat="identity", fill="gray40") +
  geom_bar(stat="identity") +
  facet_grid(.~ patient, space = 'free_x', scales = 'free_x') +  
  ylab("Mutation burden") +
  xlab("") +
  scale_y_continuous(expand=c(0,0)) +
  # annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size=6),
    # axis.text.x = element_blank(),
    axis.title.y = element_text(size=7),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size=8),
    # strip.background.x = element_blank(),
    axis.line.y = element_line(colour = "black")
  )
p_burden

# Distribution driver/immune?
gene_type_t<- table(mutClon_maf$sample_alias, mutClon_maf$gene_type)
summary(prop.table(gene_type_t,1)[,"immune"])
colSums(gene_type_t)
# driver housekeeping       immune 
# 844           20           46 
prop.table(colSums(gene_type_t))
# driver housekeeping       immune 
# 0.92747253   0.02197802   0.05054945 

# Save individual plots
#######################

ggsave("results/figs/mutClon_mutations.pdf", p_mut, width = 178, height = 265*1/3, units = "mm")
ggsave("results/figs/mutClon_mutations_all.pdf", p_mut_all, width = 178, height = 265, units = "mm")
ggsave("results/figs/mutClon_subst.pdf", p_subst, width = 178, height = 265*1/5, units = "mm")
ggsave("results/figs/mutClon_variants.pdf", p_var, width = 178/2, height = 265*1/5, units = "mm")
ggsave("results/figs/mutClon_burden.pdf", p_burden, width = 178/2, height = 265*1/5, units = "mm")
ggsave("results/figs/mutClon_vaf.pdf", p_vaf, width = 178/2, height = 265/4, units = "mm")

# Plot all in multipanel
#########################

# Process
p_mut<- p_mut + theme(
  legend.position = "none",
  axis.ticks = element_blank(),
  plot.margin = unit(c(0,0,0,0), "cm"),
  strip.text = element_blank(),
)

p_var<- p_var + 
  scale_x_continuous(expand = c(0, 0)) + 
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    plot.margin = unit(c(0,0,0,0), "cm")
  )

p_subst<- p_subst + 
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    strip.text = element_blank(),
    strip.background.x = element_blank(),
  )

p_burden<- p_burden +
  annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    strip.background.x = element_blank(),
  )

# Merge plots
p<- plot_grid(
  p_burden, NULL ,p_mut, p_var, p_subst, NULL, NULL, NULL,
  ncol = 2,
  rel_heights = c(2,3,2,4),
  rel_widths =  c(2,1),
  axis = "rlbt",
  # labels = "auto",
  # label_size = 12,
  align = "hv"
)

ggsave("results/figs/mutClon_mut_overview.pdf", p, width = 178, height = 265, units = "mm")










