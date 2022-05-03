##################################
# Filter maf & add sample alias
##################################

# Load
#######
mutClon_maf<- readRDS("temp/mutClon_maf.rds")

# Filter
########

# Exclude PM03
mutClon_maf<- mutClon_maf[mutClon_maf$patient!="PM03",]

# Exclude Sample PM02010110: outlier of 4 other eyelid samples from same patient + extremely large remnant sample!
mutClon_maf<- mutClon_maf[mutClon_maf$sample!="PM02010110",]

# Exclude HLA genes (not used in this study)
mutClon_maf<- mutClon_maf[grep("HLA.*",mutClon_maf$gene, invert = T),]

# Add sample alias
##################
samples<- sort(unique(mutClon_maf$sample))
location_short<- substr(samples,7,8)
location_short[location_short=="01"]<- "E"
location_short[location_short=="02"]<- "N"
location_short[location_short=="03"]<- "G"
location_short[location_short=="04"]<- "O"
sample_alias<- paste0(1:length(samples),location_short)
names(sample_alias)<- samples
mutClon_maf$sample_alias<- factor(sample_alias[mutClon_maf$sample], sample_alias)

# Save
######
saveRDS(mutClon_maf,"data/mutClon_maf.rds")

