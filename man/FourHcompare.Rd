\name{FourHcompare}
\alias{FourHcompare}
\title{
Tests whether multiple systems have significantly different 4H-indices
}
\description{
\code{FourHcompare} uses a PERMANOVA test on the log-ratio transformed bootstrap samples of the 4H-index from multiple different hybrid systems to determine whether the systems are significantly different in terms of the relative importance of the Union, Intersection, Gain and Loss Models.
}
\usage{
FourHcompare(boots_sets, boots_types, method='ilr', permutations = 1000)
}
\arguments{
  \item{boots_sets}{
  (Required) A matrix of bootstrapped 4H-indices from all of the systems that are to be compared. Rows are individual bootstraps, and columns are the four dimensions of the 4H-index and sigma from a particular bootstrap sample. The simplest way to generate this matrix is to use rbind() to merge the outputs of \code{FourHbootstrap} from multiple different systems.
}
\item{boots_types}{
(Required) A vector specifying which hybrid system each bootstrap sample belongs to. The order of this vector must be the same as the order of the bootstrapped samples in the matrix of bootstrapped 4H indices.
}
\item{method}{
The method used to transform the compositional 4H-index in the 4-part Aitchison Simplex to a 3-dimensional euclidean vector. Statistical tests are performed on the 3-dimensional euclidean vectors. Currently, options for transformation are 'ilr' for an isometric log-ratio transformation and 'clr' for a centered log-ratio transformation. The default is to use an isometric log-ratio transformation.
}
\item{permutations}{
The desired number of permutations to use for the PERMANOVA. The default is 1000.
}
}
\details{
\code{FourHcompare} takes bootstrapped samples of the 4H-index from multiple different hybrid systems. First, each 4H-index from each bootstrapped sample is trasnformed to a 3-dimensional euclidean vector using a log-ratio transformation. This is done using the \code{ilr} or \code{clr} functions from the \code{compositions} R package. Next, the \code{PERMANOVA} function from the \code{PERMANOVA} R package is used to test whether there are significant differences between the 4H indices of the different hybrid systems.
}
\value{
This function creates a list with contents inherited from the \code{PERMANOVA} function.
}
\references{
Van den Boogaart, K. Gerald, and Raimon Tolosana-Delgado. "“Compositions”: a unified R package to analyze compositional data." Computers & Geosciences 34.4 (2008): 320-338.

van den Boogaart, K. Gerald, et al. "Package ‘compositions’." Compositional data analysis Ver (2013): 1-40.

van den Boogaart, K. Gerald, Tolosana-Delgado, R., Bren, M. "compositions: Compositional Data Analysis" https://CRAN.R-project.org/package=compositions

Vicente-Gonzalez, L., J.L. Vicente-Villardon. "PERMANOVA: Multivariate Analysis of Variance Based on Distances and Permutations" https://CRAN.R-project.org/package=PERMANOVA

Henry, L., Wickham H., et al. "rlang: Functions for Base Types and Core R and 'Tidyverse' Features" https://CRAN.R-project.org/package=rlang

}

\examples{
#Test with enterotype dataset
library(phyloseq)
data(enterotype)
#Covert the OTU table to reads, rather than fractional abundances
otu_table(enterotype)<-round(10000*otu_table(enterotype))

#Use a subset of the enterotype dataset for one hybrid system...
#...and a subset for a second hybrid system
TS<-prune_samples(grepl('TS',sample_names(enterotype)),enterotype)
MH<-prune_samples(grepl('MH',sample_names(enterotype)),enterotype)

#Randomly assign host classes (these should be known in real hybrid microbiome datasets)
#The two parent species are assigned '1' and '3' respectively, the hybrid is assigned '2'
hybrid_statusTS<-sample(1:3,154, replace=TRUE)
hybrid_statusMH<-sample(1:3,85, replace=TRUE)

#Calculate bootstrapped samples of the 4H-index for each of the two systems
bootstrapTS<-FourHbootstrap(TS,hybrid_statusTS,0.5,5,10)
bootstrapMH<-FourHbootstrap(MH,hybrid_statusMH,0.5,5,10)

#Bind the two system 4H-index matrices together
boots_to_compare<-rbind(bootstrapTS,bootstrapMH)
system_info<-c(rep(1,5),rep(2,5))

#Perform the PERMANOVA
FourHcompare(boots_to_compare, system_info, method='ilr', permutations = 1000)

}
