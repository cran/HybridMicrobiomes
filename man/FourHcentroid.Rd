\name{FourHcentroid}
\alias{FourHcentroid}
\title{
Centroid of bootstrapped 4H-indices
}
\description{
Calculates the centroid of a set of 4H-indices from bootstrapped samples.
}
\usage{
FourHcentroid(boots)
}
\arguments{
  \item{boots}{(Required) Output matrix from the \code{FourHbootstrap} function
}
}
\details{
\code{FourHcentroid} averages over each column of the output matrix from \code{FourHbootstrap} and prints out the average value of each dimension (\eqn{\mathcal{U}}, \eqn{\mathcal{I}}, \eqn{\mathcal{G}}, \eqn{\mathcal{L}}) of the 4H-index. It also calculates and prints out the average value for \eqn{\sigma.
}
}
\value{
This function outputs a dataframe reporting the average values of each dimension of the 4H-index, and the average value of \eqn{\sigma}.
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

#Find the centroid of the bootstrap samples
FourHcentroid(boot_samples)

}
