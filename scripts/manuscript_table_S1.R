# maf file
###########
maf<- readRDS("data/mutClon_maf.rds")
maf_df<- as.data.frame(maf)

# SimSen results
################
SimSen_df<- readRDS(file="data/SimSen_df.rds")

# Add sample alias
sample_alias<- as.character(maf$sample_alias)
names(sample_alias)<- maf$sample
sample_alias<- sample_alias[!duplicated(names(sample_alias))]

# Add sample alias for samples not used in TS
sample_alias_xtr<- c("N18", "G19", "O20")
names(sample_alias_xtr)<- c("PM02010201", "PM02010301", "PM02010404")
sample_alias<- c(sample_alias, sample_alias_xtr)

# Change
SimSen_df$sample<- sample_alias[SimSen_df$sample]

# # Check
# table(SimSen_df$sample, SimSen_df$location)
# tapply(SimSen_df$vaf, SimSen_df$location, "median")

# Bed file TS + reference!
###########################

TS_gene_set_df <- readRDS("data/TS_gene_set_df.rds")

# HLA genes not used in this study, exclude
TS_gene_set_df <- TS_gene_set_df[grep("HLA",TS_gene_set_df$gene, invert = T),]

table(TS_gene_set_df$gene_type)
# driver housekeeping       immune 
# 101           20           32 

# Add census OG/TSG information
census<- as.data.frame(readr::read_csv("downloads/cosmic/CGC/v91/cancer_gene_census.csv"))
rownames(census)<- census$`Gene Symbol`
TS_gene_set_df$OG_TSG<- census[TS_gene_set_df$gene,"Role in Cancer"]
TS_gene_set_df$OG_TSG[grep("oncogene|TSG", TS_gene_set_df$OG_TSG,invert = T)]<- NA
TS_gene_set_df$OG_TSG<- gsub(", fusion","",TS_gene_set_df$OG_TSG)

# Selection
############

load("results/data/manuscript_selection.RData")
dNdS_skin<- dNdS$PM$skin[,c('n','s','nm','nn','NS','NmS','NnS','dNdS','dNmdS','dNndS','p_dNdS','p_dNmdS','p_dNndS','q_dNdS','q_dNmdS','q_dNndS')]
dNdS_oral<- dNdS$PM$oral[,c('n','s','nm','nn','NS','NmS','NnS','dNdS','dNmdS','dNndS','p_dNdS','p_dNmdS','p_dNndS','q_dNdS','q_dNmdS','q_dNndS')]
PP2_skin<- PP2$PM$skin
PP2_oral<- PP2$PM$oral
colnames(PP2_skin)<- paste0("PP2_", colnames(PP2_skin))
colnames(PP2_oral)<- paste0("PP2_", colnames(PP2_oral))

selection_skin<- merge(dNdS_skin, PP2_skin, by=0)
selection_oral<- merge(dNdS_oral, PP2_oral, by=0)
colnames(selection_skin)[1]<- "gene"
colnames(selection_oral)[1]<- "gene"

# Save
######
WriteXLS::WriteXLS(c("TS_gene_set_df", "maf_df", "SimSen_df", "selection_skin", "selection_oral"), ExcelFileName = "results/tables/Table_S1.xls", SheetNames = c("TS_geneset", "maf", "SimSen", "selection_skin", "selection_oral"), row.names = F, col.names = T)
