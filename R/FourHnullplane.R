FourHnullplane <- function(boots, col='red') {

  rlang::englue("var: {{ col }}")

  #find the average value of sigma -> the fraction of parental microbial taxa that are found on both progenitors
  sigma<-mean(boots[,5])

  #the color for plotting the points
  colorme<-col

  #Points defining one edge of the null plane (variable G, L = 0)
  baseline<-c()
  #For a grid along G from 0 to 1...
  for (i in 1:11){
    #Define the G grid
    g<-0.1*(i-1)
    #Set L as zero
    l<-0
    #Define theta for each grid point
    theta<-g+l
    #Find the value of U_null at each grid point
    u<-(1-sigma)*(1-theta)
    #Find the value of I_null at each grid point
    ins<-sigma*(1-theta)
    #make a vector of I_null, U_null, G and L
    temp<-cbind(ins,u,g,l)
    #add the grid point information to the list
    baseline<-rbind(baseline,temp)
  }

  #Points defining a second edge of the null plane (variable L, G = 0)
  baseline2<-c()
  #For a grid along L from 0 to 1...
  for (i in 1:11){
    #Define the L grid
    l<-0.1*(i-1)
    #Set G as zero
    g<-0
    #Define theta for each grid point
    theta<-g+l
    #Find the value of U_null at each grid point
    u<-(1-sigma)*(1-theta)
    #Find the value of I_null at each grid point
    ins<-sigma*(1-theta)
    #make a vector of I_null, U_null, G and L
    temp<-cbind(ins,u,g,l)
    #add the grid point information to the list
    baseline2<-rbind(baseline2,temp)
  }

  #Transform the edge at L=0 variable G to the Aitchison Simplex
  ibaseline<-compositions::ipt(baseline)
  #Transform the edge at G=0 variable L to the Aitchison Simplex
  ibaseline2<-compositions::ipt(baseline2)

  #Pull off the first point on the null plane at L = 0, G = 0
  ptA<-ibaseline[1,]
  #Pull off the second point on the null plane at L = 0, G = 1
  ptB<-ibaseline[11,]
  #Pull off the third point on the null plane at L = 1, G = 0
  ptC<-ibaseline2[11,]

  #Put the three points of the null plane into a vector
  vvv<-rbind(ptA,ptB,ptC)
  #Use the vector to plot the plane based on the three points
  rgl:: triangles3d(vvv[,1],vvv[,2],vvv[,3],color=colorme,alpha=0.2,lit=FALSE)


}
