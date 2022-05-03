# Calculate intersction surface of 2 circles
# Cfr https://diego.assencio.com/?index=8d6ca3d82151bad815f78addf9b5c1c6
calc_is_surface_circle<- function(d, r1, r2){
  if(d>(r1+r2)) A<- 0
  else if(r2>(d+r1)) A<- r1^2*pi
  else if(r1>(d+r2)) A<- r2^2*pi
  else{
    d1<- (r1^2-r2^2+d^2)/(2*d)
    d2<- d-d1
    A<- (r1^2 * acos(d1/r1)) - (d1*sqrt(r1^2-d1^2)) + (r2^2 * acos(d2/r2)) - (d2*sqrt(r2^2-d2^2))
  }
  return(A)
}
