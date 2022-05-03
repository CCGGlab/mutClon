# Load data
library(readxl)
IM<- read_excel("/home/labgroups/ccgg/downloads/pub/martincorena_2015/aaa6806-Martincorena-SM.data.set.S1.xlsx", skip = 23)

# Check samples
IM$pt<- substr(IM$sampleID,1,7)
IM$sample<- substr(IM$sampleID,8,nchar(IM$sampleID))
sample_t<- table(IM$pt,IM$sample) # 4 pts, correct
# a ab ad af ah ai aj ak al an ap ar at av ax az  b bb bc bd be bf bg bh bi bj bk bl bm bn bo bp bq br bs bt
# PD13634 30  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# PD18003  0  0  7 18  0  0 12  0  8 26 21  0 25 32  8 16 12  0 16 11 15  7  8 14 12 27 41 19 13 15 18 15  0 14 42 14
# PD20399  0 19 12  9 10  0 11  0 13  5  9 16  4 14 15 10  0 13 17 12 24 12 15  6  7 12 13 11  0 12  6  7  7 11  0  0
# PD21910  0  8  7  0  9 22  5 11  0  0  0  0  0  0  0  0  9  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# 
# bv bw bx by bz  c cb cc cd ce cf cg ch ci cj ck cl cm cn co cp cq cr cs ct cu cv cw cx cy cz  d db dc dd de
# PD13634  0  0  0  0  0 42  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 30  0  0  0  0
# PD18003 16 10 22 14 22 11 22 20 19 32 27 34 19 21  0 16 18  9  0 14 21 20 21 31  0 12 12 16 15  8  0 22 29 19 29  5
# PD20399  0 15 11 14  0  4 14  7 10 11 23 17  0 14 18 17 18 13 20 12  7  4  9 12 11 11 10  0  8 11 10  8  6 12 15 13
# PD21910  0  0  0  0  0  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
# 
# df dg dh di dj dk dl dm dn do dp dq ds  e  f  g  h  i  j  k  l  m  n  o  p  q  r  s  t  u  v  w  x  y  z
# PD13634  0  0  0  0  0  0  0  0  0  0  0  0  0 41 38 30 32 37 44 42 28 21 33 33 52 23 20 49 30 36 55 53 46  0  0
# PD18003 17 16 11 23  8 20  0 13 29 14 16 27 13 17 21  9 17 19 20 25 26  9 22 21 24 20 16 21  7 19  0 21  0 29 14
# PD20399  4  8 10  6  6  9  6  7  8  2  0  0  0  9  8  3  9  0  0  9 11  0  5 12  4 11  5  6 11 17 10  9 23 12  5
# PD21910  0  0  0  0  0  0  0  0  0  0  0  0  0 17  3  7  8  4 13  6  4  5  9  0  9  9  7  5 10  7  7  5  8  0 10

# N biopsies? As reported in table S1
rowSums(table(IM$pt,IM$sample)>0) 
# PD13634 PD18003 PD20399 PD21910 
# 24      92      90      28 

# No information found on location. Only information available is from S1: numebr of samples
biopsy_size<- c(0.79,1.57,2.36,3.14,3.93,4.71)
n_biopsies<- data.frame(
  PD13634=c(3,3,3,9,3,3),
  PD18003=c(0,68,0,24,0,0),
  PD20399=c(13,57,0,20,0,0),
  PD21910=c(12,6,0,10,0,0),row.names = biopsy_size)
colSums(n_biopsies) # ok
total_biopsy_size<- colSums(n_biopsies*as.numeric(rownames(n_biopsies)))
save(n_biopsies, biopsy_size, total_biopsy_size, file = "data/IM2015_biopsy_info.RData")

planned_biopsy_size<- pi*(5/2)^2 # 19.63mm2
pi*((2/2)^2) # 3.14 mm2
pi*((1/2)^2) # 0.78 mm2

# Fuse mutation data per pt
# Need sample size to convert VAF! Assume homogeneous: mut load correlates to size!
barplot(sort(sample_t[1,]))
pts<- colnames(n_biopsies)
biopsy_size_pred<- NULL
for(pt in pts){
  sample_t_tmp<- sort(sample_t[pt,])
  sample_t_tmp<- sample_t_tmp[sample_t_tmp!=0]
  sample_t_tmp[]<- rep(rownames(n_biopsies),n_biopsies[,pt])
  names(sample_t_tmp)<- paste0(pt,names(sample_t_tmp)) 
  biopsy_size_pred<- c(biopsy_size_pred,sample_t_tmp)
}
IM$biopsy_size_pred<- biopsy_size_pred[IM$sampleID]

# Estimate size of each clone
IM$clone_size_pred<- 2*IM$vaf*as.numeric(IM$biopsy_size_pred)
hist(IM$clone_size_pred, breaks = 50) # Very similar to fig. 4A

# Get PP2 scores
# avi<- cbind(Chromosome=IM$chr,Start_Position=IM$pos,End_Position=IM$pos_end,Reference_Allele=IM$ref_nt,Tumor_Seq_Allele2=IM$mut_nt)
# if(!dir.exists("temp/annovar/")) dir.create("temp/annovar/",recursive = T)
# write.table(avi,paste0("temp/annovar/IM2015_avinput.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
# system(paste0("/home/labgroups/ccgg/tools/annovar/table_annovar.pl temp/annovar/IM2015_avinput.txt /home/labgroups/ccgg/tools/annovar/humandb -buildver hg19 -out temp/annovar/IM2015_avoutput -remove -protocol refGene,ensGene,ljb26_all -operation g,g,f -nastring ."))
avo<- read.table(paste0("temp/annovar/IM2015_avoutput.hg19_multianno.txt"),header=TRUE,sep="\t",row.names=NULL,colClasses = "character",quote=NULL,fill=TRUE,na.strings = ".")
# Extract variant information
IM$variant<- avo$ExonicFunc.refGene
IM$PP2<- as.numeric(avo$Polyphen2_HVAR_score)
IM$PP2[IM$impact=="Synonymous"]<- 0
IM$PP2[IM$impact%in%c("Frameshift Deletion","Nonsense")]<- 1

# Get hg38 coordinates
# cat(paste0("chr",IM$chr,":",IM$pos,"-",IM$pos_end),sep="\n",file="temp/IM_liftover_input.bed")
hg38<- read.table("temp/hglft_genome_6882_a77df0.bed",sep=":")
IM$pos_hg38<- gsub("\\-.*","",hg38$V2)

# Mutation id per patient
IM$mut_id<- paste(IM$pt,IM$chr,IM$pos,sep="_")

# Clone, predicted r
IM$clone_r_pred<- sqrt(IM$clone_size_pred/pi)

saveRDS(IM, "data/IM2015.rds")
