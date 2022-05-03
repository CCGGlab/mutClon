####################
# manuscript_UV_TS
######################

# Load
######
maf<- readRDS("data/mutClon_maf.rds")

# Remove mutations/clones that occur in multiple samples
maf<- maf[!duplicated(maf$mut_id)|is.na(maf$mut_id),]

# Add 192 subst types to data
##############################
maf$subst_type_alt<- paste0(maf$REF,">",maf$ALT)
maf$subst_type_alt<- factor(maf$subst_type_alt, levels = c("C>A","C>G","C>T","T>A","T>C","T>G", "G>T","G>C","G>A","A>T","A>G","A>C"))

# Subst types?
##############

prop.table(table(maf$subst_type[maf$organ=="skin"]))
# C>A        C>G        C>T        T>A        T>C        T>G 
# 0.13162119 0.05617978 0.56019262 0.10272873 0.11556982 0.03370787 

# % Dypyrimidines (T/C[C>T])?
maf$subst_type2<- paste0(substr(maf$subst_type3,1,1),"[",substr(maf$subst_type,1,1),">",substr(maf$subst_type,3,3),"]")

# Substitution types & transcriptional strand bias
########################################################

main_st<- c("C>A","C>G","C>T","T>A","T>C","T>G")
compl_st<- c("G>T","G>C","G>A","A>T","A>G","A>C")
strand_df_ls<- list()

for(organ in c("skin","oral")){
  strand_df<- data.frame(
    subst_type=rep(main_st,4),
    Strand=rep(c("Coding (untranscribed)","Template (transcribed)","Coding (untranscribed)","Template (transcribed)"),each=length(main_st)),
    isDiP=rep(c(T,T,F,F),each=length(main_st)),
    n=c(
      # 1) diP
      # All C>T on coding strand
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("C","T")&maf$organ==organ&maf$tx_strand=="+"])[main_st] +
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("C","T")&maf$organ==organ&maf$tx_strand=="-"])[compl_st],
      # All C>T on non-coding strand
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("C","T")&maf$organ==organ&maf$tx_strand=="-"])[main_st] +
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("C","T")&maf$organ==organ&maf$tx_strand=="+"])[compl_st],
      # 2) No diP
      # All C>T on coding strand
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("G","A")&maf$organ==organ&maf$tx_strand=="+"])[main_st] +
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("G","A")&maf$organ==organ&maf$tx_strand=="-"])[compl_st],
      # All C>T on non-coding strand
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("G","A")&maf$organ==organ&maf$tx_strand=="-"])[main_st] +
      table(maf$subst_type_alt[substr(maf$subst_type2,1,1)%in%c("G","A")&maf$organ==organ&maf$tx_strand=="+"])[compl_st]
    )
  )
  
  # % T/C in upstream?
  strand_df$prop<- NA
  for(strand in c("Coding (untranscribed)", "Template (transcribed)")){
    n_temp<- tapply(strand_df$n[strand_df$Strand==strand], strand_df$subst_type[strand_df$Strand==strand], function(x) x)
    prop_temp<- sapply(n_temp,function(x) prop.table(x)[1])
    strand_df$prop[strand_df$Strand==strand&strand_df$isDiP]<- paste0(round(100*prop_temp,0),"%")
  }
  
  # Transcriptional strand assymetry, p values
  strand_df$p<- NA
  strand_df$r<- NA
  strand_df$max<- NA # For plot
  for(st in main_st){
    n_temp<- tapply(strand_df$n[strand_df$subst_type==st], strand_df$Strand[strand_df$subst_type==st], "sum")
    poi_temp<- poisson.test(n_temp)
    strand_df$p[strand_df$subst_type==st][1]<- poi_temp$p.value
    strand_df$r[strand_df$subst_type==st][1]<- poi_temp$estimate
    strand_df$max[strand_df$subst_type==st][1]<- max(n_temp)
  }
  strand_df$sign_level<- NA
  strand_df$sign_level[strand_df$p<0.05]<- "*"
  strand_df$sign_level[strand_df$p<0.01]<- "**"
  strand_df$sign_level[strand_df$p<0.001]<- "***"
  
  # Plot
  p_strand<- ggplot(strand_df, aes(x=Strand, y=n, fill=Strand)) + 
    geom_bar(stat="identity") +
    facet_grid(.~subst_type) + 
    # facet_grid(.~subst_type, switch = "x") + # Labels bottom
    geom_text(aes(x = Strand, y = n, label = prop), vjust=0, size=(0.36*6)) +
    geom_text(aes(y = max+10, label = sign_level), x=1.5) +
    ylab("Number of mutations") +
    theme(
      legend.position = c(0.8, 0.8),
      legend.text = element_text(size=7),
      legend.title = element_text(size=8),
      legend.key.size = unit(0.5, "cm"),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text = element_text(size=6),  
      axis.title = element_text(size=7),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black")
    )
  assign(paste0("p_strand_",organ),p_strand)
  strand_df_ls[[organ]]<- strand_df 
}
# p_strand_oral
# p_strand_skin

strand_df_ls[["skin"]][1:6,]
#           subst_type                 Strand isDiP   n prop            p         r max sign_level
#   1        C>A Coding (untranscribed)  TRUE  15  60% 0.0005347215 0.4385965  57        ***
#   2        C>G Coding (untranscribed)  TRUE   8  44% 1.0000000000 1.0588235  18       <NA>
#   3        C>T Coding (untranscribed)  TRUE 169  88% 0.0538193288 1.2371795 193       <NA>
#   4        T>A Coding (untranscribed)  TRUE   4  21% 0.0022275315 0.4318182  44         **
#   5        T>C Coding (untranscribed)  TRUE  19  49% 0.5559977178 1.1818182  39       <NA>
#   6        T>G Coding (untranscribed)  TRUE   4  40% 1.0000000000 0.9090909  11       <NA>

# % diPyr?
CT_diP<- strand_df_ls$skin[strand_df_ls$skin$subst_type=="C>T",]
sum(CT_diP$n[CT_diP$isDiP])/sum(CT_diP$n)
# 87%

# x more prevalent?
sum(CT_diP$n[CT_diP$Strand=="Coding (untranscribed)"])/sum(CT_diP$n[CT_diP$Strand=="Template (transcribed)"])
#  1.237179

# Show DNVs & Tx strand bias 
#################################################
maf_dnv<- maf[!is.na(maf$variant)&maf$isDNV,]
maf_dnv$subst_type_DNV<- paste0(maf_dnv$REF,">",maf_dnv$ALT)

DNV_strand_t<- table(maf_dnv$subst_type_DNV, maf_dnv$tx_strand, maf_dnv$organ)

organ<- "skin"

n_CCTT_c<- DNV_strand_t["CC>TT","+",organ] + DNV_strand_t["GG>AA","-",organ] 
n_CCTT_t<- DNV_strand_t["CC>TT","-",organ] + DNV_strand_t["GG>AA","+",organ] 
other<- sum(DNV_strand_t[,,organ])-n_CCTT_c-n_CCTT_t

DNV_strand_df<- data.frame(
  group=factor(c("CC>TT (+)", "CC>TT (-)", "Other"), levels=c("CC>TT (+)", "CC>TT (-)", "Other")), 
  n=c(n_CCTT_c, n_CCTT_t, other)
)
DNV_strand_df$prop<- DNV_strand_df$n/sum(DNV_strand_df$n)
DNV_strand_df$yposw<- cumsum(DNV_strand_df$prop)- 0.5*DNV_strand_df$prop 

# Barplot?
p_tmp<- poisson.test(DNV_strand_df$n[DNV_strand_df$group=="CC>TT (+)"],DNV_strand_df$n[DNV_strand_df$group=="CC>TT (-)"])$p.value
p_DNV<- ggplot(DNV_strand_df, aes(x=group, y=n, fill=group)) + 
  geom_bar(stat="identity") +  
  geom_text(aes(y = max(n)+2, label = format(p_tmp, digits=3)), x=1.5, size=0.36*7) +
  ylab("Number of mutations") +
  scale_fill_manual(values = c("#F8766D", "#00BFC4","grey")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,21)) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    legend.key.size = unit(0.5, "cm"),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.y = element_line(colour = "black")
  )
p_DNV

# % CC>TT?
round(100*sum(DNV_strand_df$prop[DNV_strand_df$group%in%c("CC>TT (+)")]),1)
# 82.4

# Mutational signatures
########################

# Required libs
library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(mutSignatures)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# Preprocess
maf$chr<- paste0("chr",as.character(maf$chr)) # Match seqnames wxith BSGenome
x <- filterSNV(dataSet = maf,  seq_colNames = c("REF", "ALT")) # Filter non SNV
x$sample<- x$sample_alias

# Attach context
# x <- attachContext(mutData = x,
#                    chr_colName = "chr",
#                    start_colName = "start",
#                    end_colName = "end",
#                    nucl_contextN = 3,
#                    BSGenomeDb = hg38)
# NOT WORKING, ATTACH DIRECTLY
x$context<- as.character(getSeq(hg38,x$chr,as.numeric(x$start)-1,as.numeric(x$end)+1))

# Remove mismatches
x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")

# Compute mutType
x <- attachMutType(mutData = x,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context")

# Count
MC_counts <- countMutTypes(mutTable = x,
                           mutType_colName = "mutType",
                           sample_colName = "sample")

# Obtain COSMIC signatures
cosmic_ms<- mutSignatures::getCosmicSignatures()

# Run analysis
expo <- resolveMutSignatures(mutCountData = MC_counts,
                             signFreqData = cosmic_ms)

# Retrieve exposures (results)
exp <- expo$results$count.result
# xprt <- as.data.frame(exp, transpose = TRUE) # Data frame
xprt<- matrix(NA, nrow(exp@sampleId), nrow(exp@signatureId), dimnames = list(exp@sampleId$ID, exp@signatureId$ID)) 
for(i in 1:nrow(xprt)) xprt[i,]<- exp@exposures[[i]]

# Dominant signature?
dom<- apply(xprt, 1, function(x) order(x,decreasing = T)[1])
dom[levels(x$sample_alias)]
# 1E  2N  3N  4N  5O  6O  7E  8E  9E 10E 11N 12N 13N 14N 15O 16O 17O 
# 7   7   7   7  16  16  19   7   7   9   7   7   7   7   3  19   7 

# Relative contributions
ms_rel<- t(xprt/rowSums(xprt))
ms_rel<- rbind(
  ms_rel[c("COSMIC.7","COSMIC.4"),],
  Other=colSums(ms_rel[setdiff(rownames(ms_rel),c("COSMIC.7","COSMIC.4")),])
)

# Data frame + plot
ms_rel_df<- data.frame(
  sample=rep(colnames(ms_rel),each=nrow(ms_rel)),
  ms=factor(rep(rownames(ms_rel),ncol(ms_rel)),levels=rev(c("COSMIC.7","COSMIC.4","Other"))),
  prop=as.numeric(ms_rel)
)
ms_rel_df$location<- "skin"
ms_rel_df$location[grep("O",ms_rel_df$sample)]<- "oral"
ms_rel_df$patient<- NA
ms_rel_df$patient[substr(ms_rel_df$sample,1,nchar(ms_rel_df$sample)-1)%in%1:6]<- "PM01"
ms_rel_df$patient[substr(ms_rel_df$sample,1,nchar(ms_rel_df$sample)-1)%in%7:18]<- "PM02"
ms_rel_df$sample<- factor(ms_rel_df$sample, levels = levels(maf$sample_alias))
# ms_rel_df$patient<- substr(ms_rel_df$sample,1,4)
# ms_rel_df$location<- substr(ms_rel_df$sample,7,8)
# ms_rel_df$location_short<- NA
# ms_rel_df$location_short[ms_rel_df$location=="01"]<- "E"
# ms_rel_df$location_short[ms_rel_df$location=="02"]<- "N"
# ms_rel_df$location_short[ms_rel_df$location=="03"]<- "G"
# ms_rel_df$location_short[ms_rel_df$location=="04"]<- "O"
# ms_rel_df<- ms_rel_df[order(ms_rel_df$sample),]
# ms_rel_df$sample_short<- paste0(as.numeric(ms_rel_df$sample),ms_rel_df$location_short)
# ms_rel_df$sample_short<- factor(ms_rel_df$sample_short,levels=unique(ms_rel_df$sample_short))

# Some statistics
tapply(ms_rel_df$prop, list(ms_rel_df$ms,ms_rel_df$location), "mean")
# oral       skin
# Other    0.83094744 0.53915011
# COSMIC.4 0.03659067 0.06789999
# COSMIC.7 0.13246188 0.39294990

# Plot
p_ms<- ggplot(ms_rel_df, aes(x=sample, y=prop, fill=ms)) +
  geom_bar(stat="identity") +
  scale_fill_manual(name = "Mutational Signature", values = c("grey","blue","red"),labels = c("Other SBS", "SBS4", "SBS7"),guide = guide_legend(reverse = TRUE)) +
  facet_grid(. ~ patient, scales = "free", space = "free", switch="x") +
  scale_y_continuous(expand = c(0, 0)) + 
  ylab("Proportional contribution") +
  xlab("Sample") +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size=7),
    legend.title = element_text(size=8),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.placement = 'outside',
    strip.background.x = element_blank(),
    strip.text = element_text(size=8),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p_ms

# SimSen sequencing
####################

SimSen_df<- readRDS(file="data/SimSen_df.rds")

# Some statistics
tapply(SimSen_df$vaf,list(SimSen_df$location, SimSen_df$patient),"median")
#         PM01     PM02
# gluteal       NA 0.000000
# eyelid        NA 1.159508
# nose    3.524951 3.389831
# cheek   0.000000 0.000000

tapply(SimSen_df$vaf,SimSen_df$location,"median")
# gluteal   eyelid     nose    cheek 
# 0.000000 1.159508 3.389831 0.000000 

t.test(SimSen_df$vaf[SimSen_df$location=="nose"],SimSen_df$vaf[SimSen_df$location=="eyelid"])
# p-value = 0.0090

t.test(SimSen_df$vaf[SimSen_df$patient=="PM02"&SimSen_df$location=="nose"],SimSen_df$vaf[SimSen_df$patient=="PM02"&SimSen_df$location=="eyelid"])
# p-value = 0.1275

# Plot
GU_vaf_col<- c("red","blue","purple","orange")
names(GU_vaf_col)<- c("eyelid","nose","gluteal","cheek")
SimSen_df_sel<- subset(SimSen_df,SimSen_df$patient=="PM01")
p_SS<- ggplot(SimSen_df, aes(x=location, y=vaf, fill=location)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.25) +
  geom_dotplot(binaxis='y', stackdir='center', stackratio=0.5) +
  scale_fill_manual(values=GU_vaf_col) +
  ylim(0,6) +
  # facet_grid(patient ~ .) +
  ylab("Variant allele frequency (%)") +
  xlab("") +
  theme(
    legend.position = "none",
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black")
  )
p_SS

p_SS<- p_SS +
  stat_compare_means(data=SimSen_df, aes(x=location, y=vaf, fill=location),  label = "p.format", method = "t.test", comparisons = list(c("eyelid","nose")), label.y = 5, size=2)

# Individual genes/samples
##########################

library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
pdf("results/figs/manuscript_SimSen_ind.pdf")
for(gene in unique(SimSen_df$gene)){
  chr_tmp<- unique(SimSen_df[SimSen_df$gene==gene,"chr"])
  pos_tmp<- unique(SimSen_df[SimSen_df$gene==gene,"pos"])
  start_tmp<- min(pos_tmp)-13
  end_tmp<- min(pos_tmp)+13
  ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = chr_tmp)
  sTrack <- SequenceTrack(Hsapiens)
  gtrack <- GenomeAxisTrack()
  plotTracks(list(sTrack,gtrack,ideoTrack), main=gene, chromosome = chr_tmp, from = start_tmp, to = end_tmp, add53 = TRUE)
  cat(gene, start_tmp, end_tmp, sep=" ", "\n")
}
dev.off()
# DPH3 16264984 16265010 
# RPL13A 49487421 49487447 

# Bars?
p_SS_ind<- ggplot(subset(SimSen_df, patient=="PM02"),aes(x=location, y=vaf, fill=location)) +
  geom_bar(stat="identity") +
  facet_wrap(as.character(pos)~gene) +
  scale_fill_manual(values=GU_vaf_col) +
  ylab("Variant allele frequency (%)") +
  xlab("") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,6)) +
  theme(
    axis.text = element_text(size=6),  
    axis.title = element_text(size=7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.line.y = element_line(colour = "black")
  )

# Save plots
############

library(cowplot)

p0<- plot_grid(
  p_DNV,
  p_ms,
  ncol=2,
  rel_widths = c(1,2),
  axis = "rlbt",
  align = "hv",
  labels=c("B","C")
)

p1<- plot_grid(
  p_strand_skin,
  p0,
  ncol=1,
  # axis = "rlbt",
  align = "hv",
  labels=c("A",NA)
)

p2<- plot_grid(
  NULL, p_SS,
  ncol=1,
  rel_heights = c(2,1),
  labels = c("D",NA)
)

p<- plot_grid(
  p1, p2, NULL, NULL,
  ncol=2,
  rel_widths = c(2,1),
  rel_heights =  c(1,1)
)

ggsave("results/figs/manuscript_UV_TS.pdf",  p, width = 178, height = 265, units = "mm")

ggsave("results/figs/manuscript_SimSen_ind_bars.pdf",  p_SS_ind, width = 178/4, height = 265/4, units = "mm")

# Postprocessing in inkscape (import PDF)
# Remove bottom lines panel A
# Make non-DiP 50% transparant
# Move C label up
# Add barplot & SimSen figure from raw data results
# Remove boxes around legends
