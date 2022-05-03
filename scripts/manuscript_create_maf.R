#########################################
# Generate maf file used for manuscript
#########################################

# Load mutation data
#####################
maf<- readRDS(file = "temp/mut_df.rds")

# Process
#########

# Rename colnames
colnames(maf)[colnames(maf)=="seqnames"]<- "chr"
colnames(maf)[colnames(maf)=="VF"]<- "VAF"
colnames(maf)[colnames(maf)=="location"]<- "location_short"

# Reformat ALT
maf$ALT<- sapply(maf$ALT, function(x) as.character(x))

# Reformat organ/location to full names
# Sample ids not correctly formatted for organ!
maf$location<- NA
maf$location[maf$location_short=="01"]<- "eyelid"
maf$location[maf$location_short=="02"]<- "nose"
maf$location[maf$location_short=="03"]<- "gluteal"
maf$location[maf$location_short=="04"]<- "cheek"
maf$organ<- "skin"
maf$organ[maf$location_short=="04"]<- "oral"

# Distinguish Punch Biopsy from Remnant samples
###############################################
maf$biopsy<- "P"
maf$biopsy[substr(maf$sample,9,10)=="10"]<- "R" 

# Add surface area
####################

maf$surface_cm2<- NA

# Punch biopsies: calculate surface from punch diameter (0.5 cm)
maf$surface_cm2[maf$biopsy=="P"]<- pi*0.25^2

# Remnant samples: calculate surface using imageJ on images
samples<- c('PM01010110','PM01010210','PM02010110','PM02010210')
imageJ_matrix<- as.data.frame(matrix(NA,length(samples),10, dimnames=list(samples,c("totaal","P5_1", "P5_2","P5_3","P1_1","P1_2","P1_3","P1_4","remnant_pixel","remnant_surface"))))

# "PM01010110"
ij<- read.table("raw/donors/stalen/PM01/foto's biopten/oog li/imageJ/Results.csv", sep=",", header=T)
imageJ_matrix["PM01010110",1:nrow(ij)]<- ij$Area

# PM02010110
ij<- read.table("raw/donors/stalen/PM02/PM02 bis/oog/imageJ/Results.csv", sep=",", header=T)
imageJ_matrix["PM02010110",1:nrow(ij)]<- ij$Area

# PM02010210
ij<- read.table("raw/donors/stalen/PM02/PM02 bis/neus/imageJ/Results.csv", sep=",", header=T)
imageJ_matrix["PM02010210",1:nrow(ij)]<- ij$Area

# # "PM01010210"
ij<- read.table("raw/donors/stalen/PM01/foto's biopten/neus/imageJ/Results.csv", sep=",", header=T)
imageJ_matrix["PM01010210",1:4]<- ij$Area[1:4]

# Surface
imageJ_matrix$remnant_pixel<- imageJ_matrix$totaal-rowSums(imageJ_matrix[,2:7], na.rm=T)
imageJ_matrix$remnant_surface<- 0.25*0.25*pi*(imageJ_matrix$remnant_pixel/rowMeans(imageJ_matrix[,2:4],na.rm=T))

# # "PM01010210": hard to measure, calibrate to petri dish size instead
imageJ_matrix["PM01010210","remnant_surface"]<- 2*2*pi*(imageJ_matrix["PM01010210","remnant_pixel"]/ij$Area[5])

# Add to maf
for(i in 1:nrow(imageJ_matrix)){
  maf$surface_cm2[maf$sample==rownames(imageJ_matrix)[i]]<- imageJ_matrix$remnant_surface[i]   
}

# Manual curation of longer indels and dinucleotide variants
#############################################################

# Motivation: Shearwater doesn't recognize indels, DNVs, ...; some muts are falsely labelled as single bp frameshift deletion
# Manual check & curation of SNVs closer than 10 bps together in IGV 
maf$toCheck<- F # Check manually if less than 10 bps from each other

# Select muts closer than 10 bps from each other in same sample
for(s in unique(maf$sample)){
  maf_tmp<- maf[maf$sample==s,]
  for(g in unique(as.character(maf_tmp$gene))){
    maf_tmp2<- maf_tmp[maf_tmp$gene==g,]
    start_tmp<- sort(maf_tmp2$start)
    for(i in 1:length(start_tmp)){
      d_up<- 100 # dummy values for i==1|length(start_tmp)
      d_down<- 100
      if(i>1) d_up<- abs(start_tmp[i]- start_tmp[i-1])
      if(i<length(start_tmp)) d_down<- abs(start_tmp[i]- start_tmp[i+1])
      if(d_up<=10|d_down<=10) maf$toCheck[maf$sample==s&maf$gene==g&maf$start==start_tmp[i]]<- T
    }
  }
}
# sum(maf$toCheck, na.rm=T) # 283

# # Export & curate in IGV
maf_to_check<- maf[maf$toCheck,]
maf_to_check<- maf_to_check[order(maf_to_check$sample, as.numeric(as.character(maf_to_check$chr)), maf_to_check$start, maf_to_check$end),]
maf_to_check_processed<- maf_to_check[1,]
for(i in 2:nrow(maf_to_check)){
  start_tmp<- maf_to_check$start[i]
  end_tmp<- maf_to_check$end[i-1]
  if(abs(start_tmp-end_tmp)<=1){
    i_track<- i
    next
  } 
  else{
    maf_to_check_processed<- rbind(maf_to_check_processed, maf_to_check[i,])
    if(!is.na(i_track)) maf_to_check_processed$end[nrow(maf_to_check_processed)-1]<- maf_to_check[i-1,"end"]
    i_track<- NA
  }
}
# Check: ok
# maf_to_check_gr<- makeGRangesFromDataFrame(maf[maf$toCheck,], keep.extra.columns=T)
# maf_to_check_gr <- unlist(reduce(split(maf_to_check_gr, ~sample)))
# start(maf_to_check_gr)==maf_to_check_processed$start
# end(maf_to_check_gr)==maf_to_check_processed$end

# Add sequence
maf_to_check_processed$ref_curated<- as.character(getSeq(Hsapiens, paste0("chr",maf_to_check_processed$chr), maf_to_check_processed$start,maf_to_check_processed$end))

# Save
# write.table(maf_to_check_processed[,c("sample","gene","start","end","VAF","ref_curated")], "temp/flagged_muts3.txt", sep="\t", quote=F, row.names = F)

# Import curated file
maf_cur<- readxl::read_excel("temp/flagged_muts3_curated.xlsx")
# maf_cur$start==maf_to_check_processed$start

# Merge with rest
maf_to_check_processed$REF<- maf_cur$ref_curated
maf_to_check_processed$ALT<- maf_cur$alt_curated
maf_to_check_processed$start<- maf_cur$start_curated
maf_to_check_processed$end<- maf_cur$end_curated
maf_to_check_processed$variant_curated<- maf_cur$variant_curated
maf_to_check_processed$isCurated<- maf_cur$isCurated
maf_to_check_processed<- maf_to_check_processed[maf_to_check_processed$isCurated==1,]

# Add to main maf
maf<- maf[!maf$toCheck,]
maf$ref_curated<- NA
maf$variant_curated<- NA
maf$isCurated<- 0
maf<- rbind(maf, maf_to_check_processed)

# Update width
maf$width<- 1+maf$end-maf$start
  
# Some "frameshift substitution" are in fact fs deletions e.g. CT>C
maf$start[maf$width==2&!maf$isCurated]<- maf$start[maf$width==2&!maf$isCurated] + 1
maf$REF[maf$width==2&!maf$isCurated]<- substr(maf$REF[maf$width==2&!maf$isCurated],2,2) 
maf$ALT[maf$width==2&!maf$isCurated]<- "-"

# Annotate mutations using annovar
###################################

# avi<- cbind(
#   Chromosome=paste0('chr',maf$chr),
#   Start_Position=maf$start,
#   End_Position=maf$end,
#   Reference_Allele=maf$REF,
#   Tumor_Seq_Allele2=maf$ALT
# )
# if(!dir.exists("temp/annovar/")) dir.create("temp/annovar/",recursive = T)
# write.table(avi,paste0("temp/annovar/","mutClon","_avinput.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
# system(paste0("/home/labgroups/ccgg/tools/annovar/table_annovar.pl temp/annovar/","mutClon","_avinput.txt /home/labgroups/ccgg/tools/annovar/humandb -buildver hg38 -out temp/annovar/","mutClon","_avoutput -remove -protocol refGene,ensGene,ljb26_all -operation g,g,f -nastring ."))
avo<- read.table(paste0("temp/annovar/","mutClon","_avoutput.hg38_multianno.txt"),header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")

# Add variant
maf$variant<- avo$ExonicFunc.refGene

# Add aa change
avo_aa<- avo$AAChange.refGene
avo_aa<- strsplit(unlist(lapply(avo_aa, function(x) unlist(strsplit(x,","))[1])),":") # Take first variant
maf$RefSeq<- sapply(avo_aa, function(x) x[2])
maf$nt_change<- sapply(avo_aa, function(x) x[4])
maf$aa_change<- sapply(avo_aa, function(x) x[5])

# Check genes
sum(sapply(avo_aa, function(x) x[1])!=maf$gene, na.rm=T) # 0, perfect
# NAs? unanotated, same genes in annovar

# Check variants
mut_aa<- gsub("p\\.","",maf$aa_change)
mut_aa<- gsub("[[:digit:]]+","/",mut_aa)
maf$ref_aa<- substr(mut_aa,1,1)
maf$alt_aa<- substr(mut_aa,nchar(mut_aa),nchar(mut_aa))

# Add PP2
maf$PP2<- as.numeric(avo$Polyphen2_HVAR_score)
maf$PP2[maf$variant=="synonymous SNV"]<- 0
maf$PP2[maf$variant%in%c("frameshift deletion","frameshift insertion", "frameshift substitution","stopgain")]<- 1

# Add DNV
maf$isDNV<- !is.na(maf$variant_curated)&maf$variant_curated=="DNV"

# Reannotate DNV: all doublechecked in IGV! Cfr screenshots
View(maf[!is.na(maf$variant)&maf$isDNV,])
# Oe stopgain, rest is always nonsynonymous
maf[!is.na(maf$variant)&maf$isDNV&maf$variant=="nonframeshift substitution","variant"]<- "nonsynonymous SNV"

# Only keep relevant columns
selected_variables<- c(
  'sample',
  'patient',
  'organ',
  'location',
  'gene',
  'chr',
  'start',
  'end',
  'width',
  'strand',
  'REF',
  'ALT',
  'FW',
  'BW',
  'FD',
  'BD',
  'VAF',
  'QV',
  'RefSeq',
  'nt_change',
  'aa_change',
  'variant',
  'PP2',
  'biopsy',
  'surface_cm2',
  'isCurated',
  'isDNV'
)
maf<- maf[,selected_variables]

# Factorize all genes that have been sequenced & add source information
TS_gene_df<- readRDS("data/TS_gene_set_df.rds")
rownames(TS_gene_df)<- TS_gene_df$gene
maf$gene<- factor(maf$gene, levels = TS_gene_df$gene)
maf$gene_type<- TS_gene_df[maf$gene,"gene_type"]

# Add entrez gene id
library(org.Hs.eg.db)
maf$entrez<- mapIds(org.Hs.eg.db, keys = as.character(maf$gene), column=c("ENTREZID"), keytype='SYMBOL')

# Add transcription strand
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb<- TxDb.Hsapiens.UCSC.hg38.knownGene
# detach("package:dplyr", unload=TRUE) # Necessary to run select????
tx_strand<- select(txdb, keys = unique(maf$entrez), columns=c("TXSTRAND"), keytype="GENEID")
tx_strand[tx_strand$GENEID==tx_strand$GENEID[duplicated(tx_strand$GENEID)],"TXSTRAND"]<- NA # Abiguous strand NA
tx_strand<- tx_strand[!duplicated(tx_strand$GENEID),]
rownames(tx_strand)<- tx_strand$GENEID
maf$tx_strand<- tx_strand[maf$entrez,"TXSTRAND"]

# Add substitution types
########################
source("scripts/functions/get_subst_types.R")
st_tmp<- get_subst_types(maf$chr, maf$start, maf$ALT, isHg19 = F)
maf$subst_type<- st_tmp$subst_type
maf$subst_type3<- st_tmp$subst_type3

# Add clone size: 2*VAF*biopsy size
####################################
maf$clone_size_mm2<- 100*2*maf$VAF*maf$surface_cm2 # in mm2

# Add mut id
#############
maf$mut_id<- paste(maf$patient, maf$location, maf$gene, maf$start, sep="_")  

# Save
######
saveRDS(maf,"temp/mutClon_maf.rds")
