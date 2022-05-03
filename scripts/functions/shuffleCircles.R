# function to shuffle overlapping circles
# 
# # Example
# co1<- c(1,2)
# r1<- 3
# co2<- c(5,5)
# r2<- 1
# 
# # Plot to illustrate
# MASS::eqscplot(NA,NA,xlim=c(-6,6),ylim=c(-6,6)) # Need this for correct aspect ratio
# plotrix::draw.circle(co1[1],co1[2],r1)
# plotrix::draw.circle(co2[1],co2[2],r2)
# co2_new<- shuffleCircles(co1,co2,r1,r2)
# plotrix::draw.circle(co2_new[1],co2_new[2],r2,col="red")

shuffleCircles<- function(co1,co2,r1,r2){
  # Get d (distance between midpoints)
  dx<- co2[1]-co1[1]
  dy<- co2[2]-co1[2]
  d<- sqrt(dx^2 + dy^2) # Eucledian distance?
  # Get alpha between d and vertical lines
  a<- asin(abs(dx)/d)
  # Only shuffle if overlap
  if(d>(r1+r2)) co2_new<- co2
  else{
    # Shuffle outside if midpoint circle 2 not in circle 1
    if(d>max(r1,r2)){
      # Get coordinates when d increases to r1+r2
      delta_d<- (r1+r2)-d
      delta_x<- delta_d*sin(a)*sign(dx) # Not sure for sign, pragmatic check! Should be a way o get this in angle
      delta_y<- delta_d*cos(a)*sign(dy)
      co2_new<- co2 + c(delta_x,delta_y)
    } 
    # Shuffle inside if midpoint circle 2 in circle 1
    if(d<=max(r1,r2)){
      # Get coordinates when d increases to r1+r2
      d_new<- (r1-r2)
      x_new<- d_new*sin(a)*sign(dx) # Not sure for sign, pragmatic check! Should be a way o get this in angle
      y_new<- d_new*cos(a)*sign(dy)
      co2_new<- co1 + c(x_new, y_new)
    }
  }
  # Return
  return(co2_new)
}
