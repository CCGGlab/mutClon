take_virtual_biopsy<- function(maf, punch_biopsy_r, meta_biopsy_r, show_on_plot=T){
  
  # Get coordinates of virtual biopsy
  punch_biopsy_co<- c(
    sample(seq(punch_biopsy_r,meta_biopsy_r-punch_biopsy_r,by = 0.001),1),   
    sample(seq(punch_biopsy_r,meta_biopsy_r-punch_biopsy_r,by = 0.001),1)
  )
  if(show_on_plot) plotrix::draw.circle(punch_biopsy_co[1],punch_biopsy_co[2],punch_biopsy_r,col = rgb(0.5,0.5,0.5,0.5))    
  
  # Get mutations from sample
  maf_tmp<- maf
  maf_tmp$d<- sqrt((punch_biopsy_co[1]-maf_tmp$x)^2 + (punch_biopsy_co[2]-maf_tmp$y)^2) # Eucledian distance?
  maf_tmp$clone_size_sim<- apply(maf_tmp, 1, function(x) calc_is_surface_circle(as.numeric(x["d"]),as.numeric(x["clone_r"]),punch_biopsy_r))
  maf_tmp<- maf_tmp[maf_tmp$clone_size_sim>0,]
  
  # Recalculate vaf for biopsy
  maf_tmp$vaf_sim<- (maf_tmp$clone_size_sim/(pi*(punch_biopsy_r^2)))/2
  
  # Return
  return(maf_tmp)
}

