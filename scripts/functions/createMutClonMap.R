createMutClonMap<- function(maf, geneCols=NULL, geneBorderCols=NULL, meta_biopsy_r){
  
  # Set cols
  if(is.null(geneCols)){
    geneCols<- c("#C3A0D3FF","#C3A0D3FF","#C3A0D3FF","#D0A9ADFF","#F9D0A4FF","#ED99A4FF","#B0C1D2FF")
    names(geneCols)<- c("NOTCH1","NOTCH2","NOTCH3","TP53","FGFR3","FAT1","RBM10")
  }
  geneCols<- geneCols[unique(maf$gene_name)]
  maf$col<- NA
  maf$col<- geneCols[maf$gene_name]
  
  # Set border cols
  if(is.null(geneBorderCols)){
    geneBorderCols<- c("#5C2685FF","#5C2685FF","#5C2685FF","#882E36FF","#F38928FF","#CE2125FF","#3A618FFF")
    names(geneBorderCols)<- c("NOTCH1","NOTCH2","NOTCH3","TP53","FGFR3","FAT1","RBM10")
  }
  maf$border<- "#797B7CFF"
  maf$border<- geneBorderCols[maf$gene_name]
  geneBorderCols<- geneBorderCols[unique(maf$gene_name)]
  
  # Plot
  MASS::eqscplot(NA,NA,xlim=c(0,meta_biopsy_r),ylim=c(0,meta_biopsy_r),xlab="mm", ylab="mm", frame.plot=F) # Need this for correct aspect ratio
  for(i in 1:nrow(maf)){
    plotrix::draw.circle(maf$x[i], maf$y[i], maf$clone_r[i], col = maf$col[i], border = maf$border[i])
  }
  legend("right",legend = names(geneCols), fill = geneCols, bty="n")
}

