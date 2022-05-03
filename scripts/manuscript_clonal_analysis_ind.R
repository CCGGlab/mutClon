####################################################################
# Visualize clonaes for skin a oral & PM01 & PM02 individually
####################################################################

# Load results
load("results/data/manuscript_selection.RData")
library(ggrepel)

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

# Select genes to visualize
##########################
genes_to_plot<- c("NOTCH1", "TP53", "FAT1", "NOTCH2")
geneCols<- c("#C3A0D3FF","#D0A9ADFF","#F9D0A4FF","#ED99A4FF",NA)
names(geneCols)<- c(genes_to_plot,"Other")
geneBorderCols<- c("#5C2685FF","#882E36FF","#F38928FF","#CE2125FF","grey")
names(geneBorderCols)<- c(genes_to_plot,"Other")

# Seperate for skin & oral
############################

for(org in c("skin","oral")){
  
  # Organ
  mutClon_maf_org<- mutClon_maf[mutClon_maf$organ==org,]
  
  # For 2 donors individually
  for(pt in c("PM01", "PM02")){
    mutClon_maf_tmp<- mutClon_maf_org[mutClon_maf_org$patient==pt,]

    # Get cummulative biopsy size
    sample_surf_t<- table(mutClon_maf_tmp$sample_alias, mutClon_maf_tmp$surface_cm2)>0
    total_biopsy_size_PM<- 100*sum(as.numeric(colnames(sample_surf_t)) * t(sample_surf_t)) # in mm2
    
    # Plot
    maf_co<- getMutClonMapCoordinates(maf = mutClon_maf_tmp, total_biopsy_size = total_biopsy_size_PM, avoid_partial_overlaps = T, show_plot = T)
    maf_co$gene_name<- as.character(maf_co$gene)
    maf_co$gene_name[!maf_co$gene_name%in%genes_to_plot]<- "Other"
    
    # Don't override, maf_co = random!
    pdf(paste0("results/figs/mutClon_clonal_map_",pt,"_",org,".pdf"))
    createMutClonMap(maf = maf_co, geneCols = geneCols, geneBorderCols = geneBorderCols, meta_biopsy_r = 5) # size of 5mm
    dev.off()
    
  }
  
  
}



