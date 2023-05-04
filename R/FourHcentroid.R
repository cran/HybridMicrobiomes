FourHcentroid <- function(boots) {
  centroid<-c(mean(boots[,1]),mean(boots[,2]),mean(boots[,3]),mean(boots[,4]),mean(boots[,5]))
  intersection<-centroid[1]
  union<-centroid[2]
  gain<-centroid[3]
  loss<-centroid[4]
  sigma<-centroid[5]
  centroiddf<-data.frame(intersection,union,gain,loss,sigma)
  return(centroiddf)

}
