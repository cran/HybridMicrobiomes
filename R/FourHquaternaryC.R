FourHquaternaryC <- function(boots,col='red',addplot=FALSE,plotgrid=TRUE,size_centroid=15,size_font=2) {

  #Code for resetting3D plot
  gsi.reset3D <- function(userMatrix=diag(rep(1,4))) {
    rgl::clear3d()
    rgl::par3d(userMatrix=userMatrix)
  }

  #Modified code from the 'compositions' R package to correct an error in that code for adding to an existing plot
  plot3D.acomp2 <- function(x,parts=1:min(ncol(X),4),...,lwd=2,axis.col="gray",add=FALSE,cex=2,vlabs=colnames(x),vlabs.col=axis.col,center=FALSE,scale=FALSE,log=FALSE,bbox=FALSE,axes=TRUE,size=cex,col=1) {
    requireNamespace("rgl")
    ddd = TRUE
    X<-x
    out = NULL
    if( length(parts) == 3 ) {
      if( log ) {
        x <- clr(scale(acomp(X,parts=parts),center=center,scale=scale))
        if( ! add ) {
          gsi.reset3D()
          if( axes )
            arrows3D(diag(c(0,0,0)),diag(c(1,1,1)),labs=vlabs,col=axis.col)
        }
        rgl::points3d(x[,1, drop=ddd],x[,2, drop=ddd],x[,3, drop=ddd],size=size,...,col=col)
        out = rmult(x[,1:3, drop=ddd])
      } else {
        x <- scale(acomp(X,parts=parts),center=center,scale=scale)
        if( ! add ) {
          gsi.reset3D()
          corners <- rbind(diag(rep(1,3)),c(0,0,0))
          cl <- corners[c(1,2,3,4,1,3,2,4),]
          if( axes )
            rgl::lines3d(cl[,1],cl[,2],cl[,3],col=axis.col,size=lwd)
          if( !is.null(vlabs) )
            rgl::texts3d(corners[,1],corners[,2],corners[,3],c(vlabs,"0"),col=vlabs.col)
        }
        rgl::points3d(x[,1, drop=ddd],x[,2, drop=ddd],x[,3, drop=ddd],size=size,...,col=col)
        out = rmult(x[,1:3, drop=ddd])
      }
      rgl::rgl.viewpoint(45,35.4)
    } else if( length(parts)==4 ) {
      x <- clo(X,parts=parts)
      if( log ) {
        if( ! add ) {
          gsi.reset3D()
          corners <- normalize(ilr(diag(rep(0.5,4))+0.1))
          if( axes )
            arrows3D(corners*0,corners,col=axis.col,size=lwd,labs=vlabs)
        }
        ilrx <- ilr(scale(acomp(x),center=center,scale=scale))
        rgl::points3d(ilrx[,1, drop=ddd],ilrx[,2, drop=ddd],ilrx[,3, drop=ddd],size=size,...,col=col)
        out = rmult(ilrx[,1:3, drop=ddd])
      } else {
        corners <- diag(rep(1,4))
        if( ! add ) {
          gsi.reset3D()
          corners <- diag(rep(1,4))
          cornerlines <- corners[c(1,2,3,4,1,3,2,4),]
          cl <- ipt(cornerlines)
          if( axes )
            rgl::lines3d(cl[,1, drop=ddd],cl[,2, drop=ddd],cl[,3, drop=ddd],col=axis.col,size=lwd)
        }
        iptx <- ipt(scale(acomp(x),center=center,scale=scale))
        rgl::points3d(iptx[,1, drop=ddd],iptx[,2, drop=ddd],iptx[,3, drop=ddd],size=size,...,col=col)
        out = rmult(iptx[,1:3, drop=ddd])
        if( !is.null(vlabs) ) {
          cc <- ipt(corners)
          rgl::texts3d(cc[,1, drop=ddd],cc[,2, drop=ddd],cc[,3, drop=ddd],c(vlabs),col=vlabs.col)
        }
      }
    } else
      stop("Wrong number of parts")
    if( bbox )
      rgl::bbox3d()
    invisible(out)
  }

  rlang::englue("var: {{ boots }}")
  rlang::englue("var: {{ col }}")
  rlang::englue("var: {{ plotgrid }}")
  rlang::englue("var: {{ size_centroid }}")
  rlang::englue("var: {{ size_font }}")
  rlang::englue("var: {{ addplot }}")


  #Define the colors of the points
  colorme<-col

  #Find the centroid
  CentroidH<-c(mean(boots[,1]),mean(boots[,2]),mean(boots[,3]),mean(boots[,4]))

  #the acomp() function requires at least two points... here both points are the centroid
  DoubleCentroidH<-rbind(CentroidH,CentroidH)

  #Transform the data to the Aitchison Simplex using the acomp() function
  AcompH  <- compositions::acomp(DoubleCentroidH)

  #Define the corners of the simplex
  corners <- diag(rep(1,4))
  zz<-compositions::ipt(corners)


  #If adding to an existing plot, use the modified code from the 'compositions' package (there is an error in that package - the function below fixes the error)
  if (addplot){
    plot3D.acomp2(AcompH,  cex=size_centroid, col=colorme,  add=addplot, log=FALSE, coors=T, bbox=F,  vlabs=c('','','',''),  scale=F, center=F, axis.col=1, axes=T)
  }
  #Otherwise, use the composition package code
  else{
    compositions::plot3D.acomp(AcompH,  cex=size_centroid, col=colorme,  add=addplot, log=FALSE, coors=T, bbox=F,  vlabs=c('','','',''),  scale=F, center=F, axis.col=1, axes=T)
  }
  rgl::texts3d(1.2*zz[,1],1.2*zz[,2],1.2*zz[,3],c('intersection','union','gain','loss'),family='sans',font=1,cex=size_font)

  #Draw a grid on the simplex if required
  if (plotgrid){

    #face 1
    m11<-c(0.5,0.5,0,0)
    m12<-c(0.5,0,0,0.5)
    m13<-c(0,0.5,0,0.5)
    m14<-c(0.5,0.5,0,0)
    halfcorners1<-rbind(m11,m12,m13,m14)
    hcl1<-compositions::ipt(halfcorners1)
    p11<-c(0.25,0,0,0.75)
    p12<-c(0.25,0.75,0,0)
    p13<-c(0,0.75,0,0.25)
    t11<-c(0.75,0,0,0.25)
    t12<-c(0.75,0.25,0,0)
    t13<-c(0,0.25,0,0.75)
    tcorners1<-rbind(t11,t12,t13,p11,p12,p13,t11)
    tcl1<-compositions::ipt(tcorners1)
    rgl::lines3d(hcl1[,1],hcl1[,2],hcl1[,3],col='lightgray',cex=10)
    rgl::lines3d(tcl1[,1],tcl1[,2],tcl1[,3],col='lightgray',cex=10)


    #face2
    m21<-c(0.5,0.5,0,0)
    m22<-c(0,0.5,0.5,0)
    m23<-c(0.5,0,0.5,0)
    m24<-c(0.5,0.5,0,0)
    halfcorners2<-rbind(m21,m22,m23,m24)
    hcl2<-compositions::ipt(halfcorners2)
    p21<-c(0.25,0.75,0,0)
    p22<-c(0.25,0,0.75,0)
    p23<-c(0,0.25,0.75,0)
    t21<-c(0.75,0.25,0,0)
    t22<-c(0.75,0,0.25,0)
    t23<-c(0,0.75,0.25,0)
    tcorners2<-rbind(p21,p22,p23,t21,t22,t23,p21)
    tcl2<-compositions::ipt(tcorners2)
    rgl::lines3d(hcl2[,1],hcl2[,2],hcl2[,3],col='lightgray',cex=10)
    rgl::lines3d(tcl2[,1],tcl2[,2],tcl2[,3],col='lightgray',cex=10)


    #face3
    m31<-c(0,0.5,0.5,0)
    m32<-c(0,0,0.5,0.5)
    m33<-c(0,0.5,0,0.5)
    m34<-c(0,0.5,0.5,0)
    halfcorners3<-rbind(m31,m32,m33,m34)
    hcl3<-compositions::ipt(halfcorners3)
    p31<-c(0,0.25,0.75,0)
    p32<-c(0,0.25,0,0.75)
    p33<-c(0,0,0.25,0.75)
    t31<-c(0,0.75,0.25,0)
    t32<-c(0,0.75,0,0.25)
    t33<-c(0,0,0.75,0.25)
    tcorners3<-rbind(p31,p32,p33,t31,t32,t33,p31)
    tcl3<-compositions::ipt(tcorners3)
    rgl::lines3d(hcl3[,1],hcl3[,2],hcl3[,3],col='lightgray',cex=10)
    rgl::lines3d(tcl3[,1],tcl3[,2],tcl3[,3],col='lightgray',cex=10)

    #face4
    m41<-c(0,0,0.5,0.5)
    m42<-c(0.5,0,0.5,0)
    m43<-c(0.5,0,0,0.5)
    m44<-c(0,0,0.5,0.5)
    halfcorners4<-rbind(m41,m42,m43,m44)
    hcl4<-compositions::ipt(halfcorners4)
    p41<-c(0.25,0,0.75,0)
    p42<-c(0.25,0,0,0.75)
    p43<-c(0,0,0.25,0.75)
    t41<-c(0.75,0,0.25,0)
    t42<-c(0.75,0,0,0.25)
    t43<-c(0,0,0.75,0.25)
    tcorners4<-rbind(p41,p42,p43,t41,t42,t43,p41)
    tcl4<-compositions::ipt(tcorners4)
    rgl::lines3d(hcl4[,1],hcl4[,2],hcl4[,3],col='lightgray',cex=10)
    rgl::lines3d(tcl4[,1],tcl4[,2],tcl4[,3],col='lightgray',cex=10)
  }





}

