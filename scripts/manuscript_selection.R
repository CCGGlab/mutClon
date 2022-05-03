#########################
# Selection signals?
#########################

# Load datasets
###############

TS_gene_df<- readRDS("data/TS_gene_set_df.rds")

# PM study 
maf_PM<- readRDS("data/mutClon_maf.rds")

# IM 2015
maf_IM2015<- readRDS("data/IM2015.rds")

# SCC
maf_SCC<- readRDS("data/SCC_Chang_2020_maf.rds")

# Simulations
load("data/PM_simulations.rds")

# Merge mafs
#############

maf<- maf_PM[,c("gene", "sample", "organ", "variant", "PP2", "gene_type","mut_id")]
maf$study<- "PM"

maf_IM2015$organ<- "skin"
maf_IM2015$gene<- maf_IM2015$gene_name
maf_IM2015$study<- "IM2015"
maf_IM2015<- maf_IM2015[maf_IM2015$gene%in%levels(maf_PM$gene),]
maf_IM2015$gene_type<- "driver"
maf<- rbind(maf, maf_IM2015[,colnames(maf)])

maf_SCC$sample<- maf_SCC$Sample
maf_SCC$organ<- "skin"
maf_SCC$gene<- maf_SCC$Hugo_Symbol
maf_SCC$study<- "SCC"
maf_SCC<- maf_SCC[maf_SCC$gene%in%levels(maf_PM$gene),]
maf_SCC$gene_type<- NA
for(type in c("driver","immune","housekeeping")) maf_SCC$gene_type[maf_SCC$Hugo_Symbol%in%TS_gene_df$gene[TS_gene_df$gene_type==type]]<- type
maf_SCC$mut_id<- NA # Not/less relevant for cancer, no sample IDs
maf<- rbind(maf, maf_SCC[,colnames(maf)])

# Filter on point mutations
maf<- maf[maf$variant%in%c("nonsynonymous SNV", "startloss", "stopgain", "stoploss", "synonymous SNV"  ),]

# Remove mutations/clones that occur in multiple samples
maf<- maf[!duplicated(maf$mut_id)|is.na(maf$mut_id),]

# Calculate dNdS (method Van den Eynden 2017)
#############################################

# Merge maf with gene types for simplicity
maf_tmp<- maf
maf_tmp$gene<- maf_tmp$gene_type
maf<- rbind(maf, maf_tmp)
maf_tmp$gene<- "all"
maf<- rbind(maf, maf_tmp)

# Get n or s mutations
maf_t_all<- table(maf$gene, maf$variant, maf$organ, maf$study)

dNdS<- list()
for(loc in c("skin","oral")){
  
  for(study in c("PM", "IM2015", "SCC")){
      
    # Data for study & specific location
    maf_t<- maf_t_all[,,loc,study]
    dNdS_df<- as.data.frame(matrix(NA,nrow(maf_t),6,dimnames = list(rownames(maf_t),c("n","s","NS","dNdS","p","q"))))
    
    # Add n & s to df
    dNdS_df$s<- maf_t[,"synonymous SNV"]
    dNdS_df$n<- rowSums(maf_t) - dNdS_df$s
    dNdS_df$nm<- maf_t[,"nonsynonymous SNV"]
    dNdS_df$nn<- maf_t[,"stopgain"]
    
    # Get NS
    NS<- NS_df[rownames(maf_t),"SKCM"]
    if(loc=="oral") NS<- NS_df[rownames(maf_t),"HNSC"]
    NS[is.na(NS)]<- mean(NS[!is.na(NS)]) # Approximation
    dNdS_df$NS<- NS

    NmS<- NmS_df[rownames(maf_t),"SKCM"]
    if(loc=="oral") NS<- NS_df[rownames(maf_t),"HNSC"]
    NmS[is.na(NmS)]<- mean(NmS[!is.na(NmS)]) # Approximation
    dNdS_df$NmS<- NmS

    NnS<- NnS_df[rownames(maf_t),"SKCM"]
    if(loc=="oral") NnS<- NnS_df[rownames(maf_t),"HNSC"]
    NnS[is.na(NnS)]<- mean(NnS[!is.na(NnS)]) # Approximation
    dNdS_df$NnS<- NnS
    
    # dN/dS
    dNdS_df$dNdS<- (dNdS_df$n+1)/(dNdS_df$s+1)/dNdS_df$NS
    dNdS_df$dNmdS<- (dNdS_df$nm+1)/(dNdS_df$s+1)/dNdS_df$NmS
    dNdS_df$dNndS<- (dNdS_df$nn+1)/(dNdS_df$s+1)/dNdS_df$NnS
    
    # P
    dNdS_df$p_dNdS<- apply(dNdS_df,1, function(x){
      if((x["n"]+x["s"])==0) p_tmp<- 1
      else p_tmp<- binom.test(x["n"], x["n"]+x["s"], p = x["NS"]/(x["NS"]+1), alternative = "greater")$p.value # One sided test, only test for positive selection
      p_tmp
    }) 

    dNdS_df$p_dNmdS<- apply(dNdS_df,1, function(x){
      if((x["nm"]+x["s"])==0) p_tmp<- 1
      else p_tmp<- binom.test(x["nm"], x["nm"]+x["s"], p = x["NmS"]/(x["NmS"]+1), alternative = "greater")$p.value
      p_tmp
    }) 

    dNdS_df$p_dNndS<- apply(dNdS_df,1, function(x){
      if((x["nn"]+x["s"])==0) p_tmp<- 1
      else p_tmp<- binom.test(x["nn"], x["nn"]+x["s"], p = x["NnS"]/(x["NnS"]+1), alternative = "greater")$p.value
      p_tmp
    }) 
    
    # q
    dNdS_df$q_dNdS<- p.adjust(dNdS_df$p_dNdS,"fdr")
    dNdS_df$q_dNmdS<- p.adjust(dNdS_df$p_dNmdS,"fdr")
    dNdS_df$q_dNndS<- p.adjust(dNdS_df$p_dNndS,"fdr")
    
    # Add to list
    dNdS[[study]][[loc]]<- dNdS_df
  }
  
}

# Calculate PP2
###############

genes<- levels(maf$gene)
PP2<- list()

# dPP2
for(study in c("PM", "IM2015", "SCC")){
  cat(study, " ")
  for(loc in c("skin", "oral")){
    PP2_df<- as.data.frame(matrix(NA,length(genes),3,dimnames = list(genes,c("obs","exp","p"))))
    for(gene in genes){
    
      # Observed
      PP2_obs<- maf$PP2[maf$study==study&maf$gene==gene&maf$organ==loc]
      if(sum(!is.na(PP2_obs))==0) next
      PP2_df[gene,"obs"]<- median(PP2_obs, na.rm=T)
      # Expected
      PP2_exp_tmp<- PP2_exp[[gene]][["SKCM"]]
      if(is.null(PP2_exp_tmp)) next
      if(loc=="oral") PP2_exp_tmp<- PP2_exp[[gene]][["HNSC"]]
      PP2_df[gene,"exp"]<- median(PP2_exp_tmp, na.rm=T)
      # P
      PP2_df[gene,"p"]<- wilcox.test(PP2_exp_tmp, PP2_obs, alternative = "less")$p.value  # One sided test, only test for positive selection
    }
    # q
    PP2_df$q<- p.adjust(PP2_df$p,"fdr")
    # Add to list
    PP2[[study]][[loc]]<- PP2_df
  }
}

# Focus on HI
HI_cu<- 0.7
maf$isHI<- maf$PP2>HI_cu
maf_t_all<- table(maf$gene, maf$isHI, maf$organ, maf$study)
PP2_HI<- list()
for(study in c("PM", "IM2015", "SCC")){
  cat(study, " ")
  for(loc in c("skin", "oral")){
    
    # Data for study & specific location
    maf_t<- maf_t_all[,,loc,study]
    PP2_HI_df<- as.data.frame(matrix(NA,nrow(maf_t),7,dimnames = list(rownames(maf_t),c("n_HI","n_LI","obs_prop","exp_prop", "dHIdLI","p","q"))))
    
    # Add HI & LI to df
    PP2_HI_df$n_HI<- maf_t[,"TRUE"]
    PP2_HI_df$n_LI<- maf_t[,"FALSE"]

    # Get prop
    PP2_HI_df$obs_prop<- (PP2_HI_df$n_HI + 1)/(PP2_HI_df$n_HI + PP2_HI_df$n_LI +1)
    
    # Get expected prop
    PP2_HI_df$exp_prop<- sapply(PP2_exp, function(x) mean(x$SKCM>HI_cu, na.rm=T))[rownames(PP2_HI_df)]
    if(loc=="oral") PP2_HI_df$exp_prop<- sapply(PP2_exp, function(x) mean(x$HNSC>HI_cu, na.rm=T))[rownames(PP2_HI_df)]
    
    # dHI/dLI
    PP2_HI_df$dHIdLI<- PP2_HI_df$obs_prop/PP2_HI_df$exp_prop

    # P
    PP2_HI_df$p<- apply(PP2_HI_df,1, function(x){
      if((x["n_HI"]+x["n_LI"])==0) p_tmp<- 1
      else if(is.na(x["exp_prop"])) p_tmp<- NA
      else p_tmp<- binom.test(x["n_HI"], x["n_HI"]+x["n_LI"], p = x["exp_prop"], alternative="greater")$p.value  # One sided test, only test for positive selection
      p_tmp
      p_tmp
    }) 
    
    # q
    PP2_HI_df$q<- p.adjust(PP2_HI_df$p,"fdr")

    # Add to list
    PP2_HI[[study]][[loc]]<- PP2_HI_df
  }
}
# View(PP2_HI$PM$skin)
# View(PP2_HI$PM$oral)

# Save
######
save(dNdS, PP2, PP2_HI, maf, file="results/data/manuscript_selection.RData")

