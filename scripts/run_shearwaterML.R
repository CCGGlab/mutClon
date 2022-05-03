# Load libraries
#################
library(deepSNV)
library(GenomicRanges)
library(readxl)
library(biomaRt)

# Get ensembl genes & Granges
###############################
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get genes used for TS
genes<- readRDS("data/TS_gene_set_df.rds")$gene
gene_ranges <- getBM(attributes = c("external_gene_name","chromosome_name",'start_position','end_position'),
                     filters = "external_gene_name",
                     values = genes,
                     mart = ensembl)

# Only canonical chromosomes
gene_ranges <- gene_ranges[gene_ranges$chromosome_name%in%c(1:22,"X","Y"),]

# Get metadata
###############
bamloc<- "temp/bam"
bam_files<- list.files(bamloc, recursive = T)
bam_files<- bam_files[grep("markDup_sorted\\.bam", bam_files)]
bam_files<- bam_files[grep("\\.bai", bam_files, invert = T)]
samples<- gsub("\\/.*","",bam_files)
patients<- substr(samples,1,4)
organ<- substr(samples,5,6)
location<-  substr(samples,7,8)

meta_df<- data.frame(
  sample=samples,
  patient=patients,
  organ,
  location,
  bam_file=paste0(bamloc,"/",bam_files)
  )

## Run ShearwaterML ##
######################

mut_df<- NULL
for (i in 1:nrow(gene_ranges)){
  cat(i, " ")
  gene <- gene_ranges$external_gene_name[i]
  regions <- GRanges(seqnames = gene_ranges$chromosome_name[i], IRanges(start = gene_ranges$start_position[i], end=gene_ranges$end_position[i]))
  # print(paste0('calculating for ', gene))
  for (patient in unique(patients)){
    cat(patient, " ")
    samples_pt<- meta_df$sample[meta_df$patient==patient]
    for(sample in samples_pt){
      cat(sample," ")
      sample_loc<- meta_df$location[meta_df$sample==sample]
      samples_exl<- setdiff(meta_df$sample[meta_df$location==sample_loc], sample)
      ControlFiles<- meta_df$bam_file[meta_df$patient==patient&!meta_df$sample%in%samples_exl]   ## use other samples from same patient as genomic control
      counts <- loadAllData(ControlFiles, regions, q=30) #mq = 30 
      
      ## ShearwaterML: “betabinLRT” calculates p-values for each possible mutation
      pvals <- betabinLRT(counts, rho=1e-4, truncate=0.005)$pvals
      qvals <- p.adjust(pvals, method="BH")
      dim(qvals) = dim(pvals)
      vcfML = qvals2Vcf(qvals, counts, regions, samples = ControlFiles, mvcf = TRUE)
      
      writeVcf(vcfML, paste0('temp/vcf/',sample, '_',gene,'.vcf'), index = FALSE)
      f<- rownames(colData(vcfML))[grep(sample, rownames(colData(vcfML)))]
      # for(f in rownames(colData(vcfML))){
        indices <- geno(vcfML)$VF[,f]>0 & geno(vcfML)$QV[,f]<0.05
        if(!is.na(sum(indices)) & sum(indices) > 0){
          gt_tmp<- sapply(geno(vcfML), function(x) x[indices,f], simplify=F)
          sample_df<- data.frame(
            sample = meta_df$sample[meta_df$bam_file==f],
            patient = meta_df$patient[meta_df$bam_file==f],
            organ = meta_df$organ[meta_df$bam_file==f],
            location = meta_df$location[meta_df$bam_file==f],
            gene,
            bam = f,
            as.data.frame((rowRanges(vcfML)))[indices,],
            GQ = gt_tmp$GQ,
            QV = gt_tmp$QV,
            VF = gt_tmp$VF,
            FW = gt_tmp$FW,
            BW = gt_tmp$BW,
            FD = gt_tmp$FD,
            BD = gt_tmp$BD
          ) 
          mut_df <- rbind(mut_df, sample_df)
        # }
    }
    }
  }
}

saveRDS(mut_df, file = "temp/mut_df.rds")
