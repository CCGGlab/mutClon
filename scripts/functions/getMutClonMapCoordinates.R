# Function that creates a map of somatic mutation clones

# maf should contain following columns 
# 1) "vaf", 
# 2) "clone_size", 
# 3) "mut_id" (pt, chr, pos)
# 
# maf<- readRDS("temp/IM2015.rds")
# maf<- maf[maf$pt=="PD18003",]
# maf$clone_size<- maf$clone_size_pred
# 
# load("temp/IM2015_biopsy_info.RData")
# total_biopsy_size<- total_biopsy_size["PD18003"]
# 
source("scripts/functions/calc_is_surface_circle.R")
source("scripts/functions/shuffleCircles.R")

getMutClonMapCoordinates<- function(maf, total_biopsy_size, avoid_partial_overlaps=F, show_plot=T){
  
  # Pool similar muts in different biopsies from same pt
  size_recalc<- tapply(maf$clone_size,maf$mut_id,"sum")
  maf$clone_size<- size_recalc[maf$mut_id]
  maf$clone_r<- sqrt(maf$clone_size/pi)
  maf<- maf[!duplicated(maf$mut_id),]
  
  # Create square metasample from all individual samples
  meta_biopsy_size<- sum(total_biopsy_size)
  meta_biopsy_r<-  sqrt(meta_biopsy_size)
  
  # Randomly assign mutations/clones to metasample
  maf<- maf[order(maf$clone_size,decreasing = T),]
  maf$x<- NA
  maf$y<- NA
  if(show_plot) MASS::eqscplot(NA,NA,xlim=c(0,meta_biopsy_r),ylim=c(0,meta_biopsy_r)) # Need this for correct aspect ratio
  for(i in 1:nrow(maf)){
    # cat(i, " ")
    # Select random coordinates
    co_tmp<- c(
      sample(seq(0,meta_biopsy_r,0.0001),1,replace = T),
      sample(seq(0,meta_biopsy_r,0.0001),1,replace = T)
    )
    r_tmp<- maf$clone_r[i]
    # Check whether clones overlap with earlier clones: becomes very slow for high MRs!!!
    if(avoid_partial_overlaps){
      if(i>1){
        c<- 1
        attempts<- 1
        major_attempts<- 1
        while(c < i){
          # cat(c," ")
          co_new<- shuffleCircles(as.numeric(maf[c,c("x","y")]),co_tmp,maf$clone_r[c],r_tmp)
          if(!identical(co_tmp,co_new)){
            co_tmp<- co_new
            attempts<- attempts + 1
            c<- 1
            if(attempts%%10==0){ # Prevent endless loops
              c<- 1
              co_tmp<- c(
                sample(seq(0,meta_biopsy_r,0.0001),1,replace = T),
                sample(seq(0,meta_biopsy_r,0.0001),1,replace = T)
              )
              major_attempts<- major_attempts + 1 
              if(major_attempts==100) c<- i # Give up finding overlaps
            }
          }
          else{
            c<- c+1
          }
        }
      }
    }
    maf$x[i]<- co_tmp[1]
    maf$y[i]<- co_tmp[2]
    if(show_plot) plotrix::draw.circle(co_tmp[1], co_tmp[2], r_tmp)
  }
  # Return
  return(maf)
}
