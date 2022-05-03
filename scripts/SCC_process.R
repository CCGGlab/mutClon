################
# SCC_process
################

# Dataset from a recent comprehensive SCC analysis on biorxiv
maf_SCC<- as.data.frame(readxl::read_excel("downloads/pub/chang_2020/media-4.xlsx", skip = 1))

# Annovar
avi<- cbind(Chromosome=maf_SCC$Chromosome,Start_Position=maf_SCC$Start_Position,End_Position=maf_SCC$End_Position,Reference_Allele=maf_SCC$Reference_Allele,Tumor_Seq_Allele2=maf_SCC$Tumor_Seq_Allele2)
write.table(avi,paste0("temp/annovar/SCC_avinput.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
system(paste0("/home/labgroups/ccgg/tools/annovar/table_annovar.pl temp/annovar/SCC_avinput.txt /home/labgroups/ccgg/tools/annovar/humandb -buildver hg19 -out temp/annovar/SCC_avoutput -remove -protocol refGene,ensGene,ljb26_all -operation g,g,f -nastring ."))
avo<- read.table(paste0("temp/annovar/SCC_avoutput.hg19_multianno.txt"),header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")

# Extract variant information
maf_SCC$variant<- avo$ExonicFunc.refGene

# PP2
maf_SCC$PP2<- as.numeric(avo$Polyphen2_HVAR_score)
maf_SCC$PP2[maf_SCC$Variant_Classification=="Silent"]<- 0
maf_SCC$PP2[maf_SCC$Variant_Classification%in%c("Frame_Shift_Del","Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation")]<- 1

# save
saveRDS(maf_SCC, "data/SCC_Chang_2020_maf.rds")
