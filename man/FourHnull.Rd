\name{FourHnull}
\alias{FourHnull}
\title{
Bootstrapped null models of the 4H-index
}
\description{
Calculates 4H-indices for bootstrapped samples of hybrid null models drawn from an object of class phyloseq containing host-associated microbial community data from both hybrid hosts and hosts of each parent species.
}
\usage{
FourHnull(x, class_grouping, core_fraction, boot_no,
              sample_no, null_model=1,replace_hosts=FALSE, reads=NULL,
              rarefy_each_step=TRUE,seed=NULL, use_microViz='yes')
}
\arguments{
  \item{x}{
  (Required) Host-associated microbial community data . This must be in the form of a phyloseq object and must contain an OTU abundance table with information on the host-associated microbial communities from both hybrid hosts and hosts of each progenitor taxon. The OTU abundance table should contain read counts, not fractional abundances.
}
\item{class_grouping}{
(Required) A numerical vector identifying the host class of each host from the OTU table. This vector must be ordered according to the OTU abundance table, and must use the following notation: Progenitor One = 1, Hybrids = 2, Progenitor Two = 3.
}
\item{core_fraction}{
(Required) The fraction of hosts that a microbial taxon must be found on to be considered part of a host class's 'core' microbiome
}
\item{boot_no}{
(Required) The desired number of bootstrap samples
}
\item{sample_no}{
(Required) The number of hosts of each host class that should be selected to generate each bootstrap sample. This number should be less than the minimum number of hosts in any given class.
}
\item{null_model}{
The particular null model to use (see below for a description of each null model). The default is 1 (averaging microbial taxon abundances from a randomly drawn individual of each progenitor species)
}
\item{replace_hosts}{
Whether or not to replace hosts during bootstrapping. The default is FALSE (i.e., all hosts in a given bootstrap are unique).
}
\item{reads}{
The number of reads to rarefy to for each microbiome sample. The default is to use the number of reads in the sample with the lowest number of reads. Any user-specified value must be lower than the number of reads in the sample with the lowest number of reads.
}
\item{rarefy_each_step}{
Whether or not to rarefy each bootstrap sample separately. The default is TRUE (i.e., perform a new rarefaction on the OTU table for each bootstrap sample). However, for cases with large numbers of samples, large numbers of reads per sample or high microbial diversity, this can be slow; thus, there is the option to rarefy the microbiome samples a single time prior to bootstrapping.
}
\item{seed}{
Optional seed for random number generator. The default is NULL (i.e., seed picked at random at time of simulation).
}
\item{use_microViz}{
Whether or not to use the microViz package to check your phyloseq object. If installation of the microViz package is problematic, this test can be skipped, but the user should verify themselves that a phyloseq object exists AND includes sample data (the sample data placeholder is necessary for the FourHnull code to work).
}
}
\details{
\code{FourHnull} generates bootstrapped samples of hybrid null models of the 4H index. As in \code{FourHbootstrap}, each iteration begins by rarefying the OTU table such that all samples have the same library size. Once the OTU table has been rarefied, a subset, \code{sample_no}, of hosts is selected from each progenitor class. Hosts can be sampled with or without replacement. The default is to sample hosts without replacement, such that all progenitors in any given bootstrap sample are unique. Next,\code{sample_no} hybrid individuals are generated according to the chosen null model. Null models are as follows:

\code{null_model= 1}: each hybrid individual is created by randomly selecting one individual from the progenitor 1 class and one individual from the progenitor 2 class and then averaging the microbial abundances from the two progenitor individuals.

\code{null_model= 2}: each hybrid individual is created by randomly selecting two individuals from the progenitor 1 class and one individual from the progenitor 2 class and then averaging the microbial abundances from the three progenitor individuals.

\code{null_model= 3}: each hybrid individual is created by randomly selecting one individual from the progenitor 1 class and two individuals from the progenitor 2 class and then averaging the microbial abundances from the three progenitor individuals.

\code{null_model= 4}: each hybrid individual is created by randomly selecting two individual from the progenitor 1 class and then averaging the microbial abundances from these two progenitor individuals.

\code{null_model= 5}: each hybrid individual is created by randomly selecting two individual from the progenitor 2 class and then averaging the microbial abundances from these two progenitor individuals.

\code{null_model= 6}: each hybrid individual is created by randomly selecting one individual from the progenitor 1 class and using its microbial abundances.

\code{null_model= 7}: each hybrid individual is created by randomly selecting one individual from the progenitor 2 class and using its microbial abundances.

\code{null_model= 10}: each hybrid individual is created by randomly selecting one individual from the progenitor 1 class and one individual from the progenitor 2 class, randomly setting half of the microbial taxon abundances only found on the progenitor 1 individual to zero, randomly setting half of the microbial taxon abundances only found on the progenitor 2 individual to zero, and then setting all other microbial taxon abundances to one.

Once null model hybrid iterates have been generated, \code{FourHnull} proceeds as in \code{FourHbootstrap} and calculates the 4H-index and \eqn{\sigma} for each bootstrap sample.

}
\value{
This function returns a matrix, where rows are individual bootstrap samples and columns are the four dimensions of the 4H-index, as well as \eqn{\sigma}.
}
\references{
Camper, B., Laughlin, Z. et al. "A Conceptual Framework for Host-Associated Microbiomes of Hybrid Organisms" <doi:10.1101/2023.05.01.538925>

Henry, L., Wickham H., et al. "rlang: Functions for Base Types and Core R and 'Tidyverse' Features" https://CRAN.R-project.org/package=rlang

McMurdie, Paul J., and Susan Holmes. "phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data." PloS one 8.4 (2013): e61217.

McMurdie, Paul J., and Susan Holmes. "Phyloseq: a bioconductor package for handling and analysis of high-throughput phylogenetic sequence data." Biocomputing 2012. 2012. 235-246.
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

FourHnull(enterotype,hybrid_status,0.5,5,10, use_microViz='no')

}
