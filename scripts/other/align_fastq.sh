#!/bin/bash

# Reference genome
refgenome_filename="downloads/genomes/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Get all fastqfiles
fastq_files=$(find raw/BGI/ | grep "\.fq" -)

# Get all samples to be processed
selected_samples=$(find raw/BGI/ | grep "L01_.*_1\." - )
selected_samples=$(cut -d'/' -f5 <<<$selected_samples) # 5th element = sample directory

# Process by sample
for sample in $selected_samples; 
# sample=($selected_samples[1])
do 

  bamdir=temp/bam/$sample

  # Check whether processed already
  if [ -d  $bamdir ];
  then
      echo "Sample processed already, skipping ...";
      continue;
  fi;

  # Create patient directory (-p if not exists)
  mkdir -p $bamdir
  
  for i in 1 2
  do
  
  # Get both fq files
  file1=$(find $fastq_files | grep $sample - | grep "L0"$i - | grep _1\. -)
  file2=$(find $fastq_files | grep $sample - | grep "L0"$i - | grep _2\. -)
  
  # Align to grch38 using bwa-mem2: only takes a few minutes
  # -t number of threads
  # DeepSNV loadAllData crashes when q parameter is not used here!
  # bwa-mem2 mem -t 24 $refgenome_filename $file1 $file2 | samtools view -@ 24 -bS - > $bamdir/"L0"$i.bam
  bwa-mem2 mem -t 24 $refgenome_filename $file1 $file2 | samtools view -@ 24 -b -S -q 20 - > $bamdir/"L0"$i.bam
  
  # RG (+sort) - downstream errors avoided
  eval 'file_tmp=("${file'$i'}")'
  file_tmp=$(echo $file_tmp | sed 's/..*\///g' -)
  # PU=$(echo $(cut -d'_' -f3 <<<$file_tmp).$(cut -d'_' -f2 <<<$file_tmp))
  PU=$(echo $(cut -d'_' -f3 <<<$file_tmp)) # Don't add lane to remove dupl afterwards!
  LB=$(cut -d'_' -f1 <<<$file_tmp)
  picard AddOrReplaceReadGroups -I $bamdir/"L0"$i.bam -O $bamdir/"L0"$i'_rg.bam' -LB $LB -PL ILLUMINA -PU $PU -SM $sample -SO coordinate -CREATE_INDEX true

  done

# Merge 
samtools merge -@ 24 $bamdir/merged.bam $bamdir/L01_rg.bam $bamdir/L02_rg.bam  

# Sort & index
samtools sort -@ 24 $bamdir/merged.bam -o $bamdir/merged_sorted.bam
# samtools index $bamdir/merged_sorted.bam

# Mark duplicates
picard MarkDuplicates -I $bamdir/merged_sorted.bam  -M $bamdir/merged_bam_metrix.txt -O $bamdir/merged_markDup.bam --REMOVE_DUPLICATES --REMOVE_SEQUENCING_DUPLICATES

# Sort & index
samtools sort -@ 24 $bamdir/merged_markDup.bam -o $bamdir/merged_markDup_sorted.bam
samtools index -@ 24 $bamdir/merged_markDup_sorted.bam

# Get bam stats: includes insert sizes, read lengths, ...
samtools stats -@ 24 $bamdir/merged_markDup_sorted.bam > $bamdir/merged_sorted_bam_stat.txt

done  
