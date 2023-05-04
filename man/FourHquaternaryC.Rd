\name{FourHquaternaryC}
\alias{FourHquaternaryC}
\title{
Plot centroid of bootstrapped 4H-indices
}
\description{
Draws the centroid of bootstrapped values of the 4H-index on a quaternary (four dimensional barycentric, Aitchison Simplex) plot.
}
\usage{
FourHquaternaryC(boots,col='red',addplot=FALSE,plotgrid=TRUE,
                                            size_centroid=15,size_font=2)
}
\arguments{
  \item{boots}{Output from the \code{FourHbootstrap} function
}
\item{col}{The color of the points in the quaternary plot. The default is 'red'}
\item{addplot}{Whether or not to add the data to an existing plot. The default is FALSE}
\item{plotgrid}{Whether or not to plot gridlines on the quaternary plot. The default is TRUE}
\item{size_centroid}{The size of the point representing the centroid of the bootstrap samples. The default is 15.}
\item{size_font}{The size of the font used to label the four verticies. The default is 2.}
}
\details{
\code{FourHquaternaryC} allows for visualization of the centroid of bootstrapped 4H-indices within the framework of a 4-dimensional Aitchison Simplex, also known as a 4-dimensional barycentric plot or quaternary plot. Commonly used for compositional data, the Aitchison Simplex is convenient for the 4H-index, because \eqn{\mathcal{U} + \mathcal{I} + \mathcal{G} + \mathcal{L} = 1}, reflecting the fact that each of the four different conceptual models - Union, Intersection, Gain and Loss - comprise one part of the whole, where the whole is the total microbial diversity across the entire hybrid host complex. \code{FourHquaternaryC} uses functions from the \code{compositions} R package to transform the four coordinates of the 4H-index, such that they can be drawn on a quaternary plot.
}
\value{
This function plots the centroid of bootstrapped 4H-indices on quaternary plot (Aitchison Simplex, 4-dimensional barycentric plot).

}
\references{
Henry, L., Wickham H., et al. "rlang: Functions for Base Types and Core R and 'Tidyverse' Features" https://CRAN.R-project.org/package=rlang

Van den Boogaart, K. Gerald, and Raimon Tolosana-Delgado. "“Compositions”: a unified R package to analyze compositional data." Computers & Geosciences 34.4 (2008): 320-338.

van den Boogaart, K. Gerald, et al. "Package ‘compositions’." Compositional data analysis Ver (2013): 1-40.

van den Boogaart, K. Gerald, Tolosana-Delgado, R., Bren, M. "compositions: Compositional Data Analysis" https://CRAN.R-project.org/package=compositions

Adler, Daniel, Oleg Nenadic, and Walter Zucchini. "Rgl: A r-library for 3d visualization with opengl." Proceedings of the 35th symposium of the interface: computing science and statistics, Salt Lake City. Vol. 35. 2003.

Adler, D., Murdoch, M. D., Suggests, M. A. S. S., WebGL, P. L. Y., OBJ, S., & OpenGL, S. (2019). Package ‘rgl’.

Murdoch, Duncan, Daniel Adler, and Oleg Nenadic. "Package ‘rgl’." R Package (2023).


}

\examples{
#Test with enterotype dataset
library(phyloseq)
data(enterotype)
#Covert the OTU table to reads, rather than fractional abundances
otu_table(enterotype)<-round(10000*otu_table(enterotype))

#Randomly assign host classes (these should be known in a real hybrid microbiome dataset)
#The two parent species are assigned '1' and '3' respectively, the hybrid is assigned '2'
hybrid_status<-sample(1:3,280, replace=TRUE)

#Bootstrap the dataset
boot_samples<-FourHbootstrap(enterotype,hybrid_status,0.5,5,10)

#Plot the centroid of the bootstrapped samples as a new quaternary plot
FourHquaternaryC(boot_samples,col='red')

#Plot the centroid of the bootstrapped samples on an existing quaternary plot
#An existing quaternary plot should already be open
FourHquaternaryC(boot_samples,col='red',addplot=TRUE)

}
