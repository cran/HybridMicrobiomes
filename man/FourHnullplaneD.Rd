\name{FourHnullplaneD}
\alias{FourHnullplaneD}
\title{
Bootstrap Distribution Around the Null Plane
}
\description{
Calculates the distribution of bootstrap samples around the null expectation for Union vs. Intersection models.
}
\usage{
FourHnullplaneD(boots)
}
\arguments{
  \item{boots}{(Required) Output matrix from the \code{FourHbootstrap} function
}


}
\details{
For each bootstrap sample, \code{FourHnullplaneD} uses the fraction of microbes shared by both progenitors to calculate null model expectations for the value of \eqn{\mathcal{U}} and \eqn{\mathcal{I}}. It then calculates the difference between the observed and expected values along the intersection and union dimensions \eqn{\mathcal{I}-\mathcal{I}_{null}} and \eqn{\mathcal{U}-\mathcal{U}_{null}}. Finally, \code{FourHnullplaneD} plots a histogram of the difference between observed and expected intersection dimensions, \eqn{\mathcal{I}-\mathcal{I}_{null}} and calculates the average distance between the observed and expected intersection dimensions, as well as the standard deviation of \eqn{\mathcal{I}-\mathcal{I}_{null}} and the fraction of bootstrapped samples where \eqn{\mathcal{I}-\mathcal{I}_{null}<0} (i.e., hybrids preferentially retain microbial taxa only found on one progenitor). Preferential retention of microbial taxa found on both progenitors can then be assessed as >95\% of bootstrapped samples having \eqn{\mathcal{I}-\mathcal{I}_{null}>0}.
}
\value{
This function outputs a list including a dataframe with the following lists:
\item{\code{sigma}}{The inputted values of \eqn{\sigma} for each bootstrapped sample.}
\item{\code{theta}}{The calculated values of \eqn{\theta = \mathcal{G}+\mathcal{L}} for each bootstrapped sample.}
\item{\code{union}}{The observed values of the union dimension, \eqn{\mathcal{U}}, for each bootstrapped sample.}
\item{\code{union_expectation}}{The expected values of the union dimension, \eqn{\mathcal{U}_{null}}, for each bootstrapped sample.}
\item{\code{intersection}}{The observed values of the intersection dimension, \eqn{\mathcal{I}}, for each bootstrapped sample.}
\item{\code{intersection_expectation}}{The expected values of the intersection dimension, \eqn{\mathcal{I}_{null}}, for each bootstrapped sample.}
\item{diffI}{The difference between the observed and expected values of the intersection dimension for each bootstrapped sample.}
\item{diffU}{The difference between the observed and expected values of the union dimension for each bootstrapped sample.}
as well as the average distance and standard deviation between the observed and expected intersection dimensions, and the fraction of bootstrapped samples where hybrids preferentially retain microbial taxa only found on one progenitor.
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

#Calculate the expected and observed distributions of bootstrapped samples around the null plane
FourHnullplaneD(boot_samples)

}
