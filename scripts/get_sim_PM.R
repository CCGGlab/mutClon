# Get expected NS ratio for genes

# Load simulated mutation data
GPPM<- readRDS("data/GPPM_PM.rds")

# Get triNT mutation probability data from TCGA
ms_all<- readRDS("data/TCGA_mutProb.rds")

################################
# Calculate expected N/S ratios
################################

# Create variable table
GPPM$var<- NA
GPPM$var[GPPM$variant%in%c("nonsynonymous SNV", "stopgain", "stoploss")]<- "n"
GPPM$var[GPPM$variant%in%c("synonymous SNV")]<- "s"
var_t<- table(GPPM$gene, GPPM$var, GPPM$subst_type3)
variant_t<- table(GPPM$gene, GPPM$variant, GPPM$subst_type3)

# Add groups: driver, immune, housekeeping genes
TS_gene_df<- readRDS("data/TS_gene_set_df.rds")
rownames(TS_gene_df)<- TS_gene_df$gene
GPPM$gene_type<- NA
for(type in c("driver","immune","housekeeping")) GPPM$gene_type[GPPM$gene%in%TS_gene_df$gene[TS_gene_df$gene_type==type]]<- type
var_type_t<- table(GPPM$gene_type, GPPM$var, GPPM$subst_type3)
variant_type_t<- table(GPPM$gene_type, GPPM$variant, GPPM$subst_type3)

# Calculate per gene & MS
genes<- c(sort(unique(rownames(var_t))),c("driver","immune","housekeeping","all"))
NS_df=NmS_df=NnS_df<- data.frame(row.names=genes)

for(ms in colnames(ms_all)){
  cat(ms," ")
  NS_df[ms]<- NA
  for(gene in genes){
    # NS
    if(gene%in%c("driver","immune","housekeeping")) var_gene_t<- var_type_t[gene,,]
    else if(gene=="all") var_gene_t<- table(GPPM$var, GPPM$subst_type3)
    else var_gene_t<- var_t[gene,,]
    N<- sum(var_gene_t["n",rownames(ms_all)]*ms_all[,ms])
    S<- sum(var_gene_t["s",rownames(ms_all)]*ms_all[,ms])
    NS<- N/S
    NS_df[gene,ms]<- NS
    # NmS & NnS
    if(gene%in%c("driver","immune","housekeeping")) variant_gene_t<- variant_type_t[gene,,]
    else if(gene=="all") var_gene_t<- table(GPPM$variant, GPPM$subst_type3)
    else variant_gene_t<- variant_t[gene,,]
    Nm<- sum(variant_gene_t["nonsynonymous SNV",rownames(ms_all)]*ms_all[,ms])
    Nn<- sum(variant_gene_t["stopgain",rownames(ms_all)]*ms_all[,ms])
    S<- sum(variant_gene_t["synonymous SNV",rownames(ms_all)]*ms_all[,ms])
    NmS<- Nm/S
    NnS<- Nn/S
    NmS_df[gene,ms]<- NmS
    NnS_df[gene,ms]<- NnS
  }  
}

#########################
# Calculate Expected PP2
##########################

PP2_exp<- list()

for(gene in genes){
  cat(gene," ")
  PP2_exp[[gene]]<- list()
  if(gene%in%c("driver","immune","housekeeping")) GPPM_gene<- GPPM[GPPM$gene_type==gene]
  else if(gene=="all") GPPM_gene<- GPPM
  else GPPM_gene<- GPPM[GPPM$gene==gene]
  for(ms in colnames(ms_all)){
    PP2_tmp<- sample(GPPM_gene$PP2,100000,T,prob = ms_all[GPPM_gene$subst_type3,ms])
    PP2_exp[[gene]][[ms]]<- PP2_tmp
  }  
}


# Save
######
save(PP2_exp, NS_df, NmS_df, NnS_df, file="data/PM_simulations.rds")


