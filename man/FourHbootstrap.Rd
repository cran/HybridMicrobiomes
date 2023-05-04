\name{FourHbootstrap}
\alias{FourHbootstrap}
\title{
Bootstraps of the 4H-index
}
\description{
Calculates 4H-indices for bootstrapped samples drawn from an object of class phyloseq containing host-associated microbial community data from both hybrid hosts and hosts of each parent species.
}
\usage{
FourHbootstrap(x, class_grouping, core_fraction, boot_no,
              sample_no, replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE)
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
\item{replace_hosts}{
Whether or not to replace hosts during bootstrapping. The default is FALSE (i.e., all hosts in a given bootstrap are unique).
}
\item{reads}{
The number of reads to rarefy to for each microbiome sample. The default is to use the number of reads in the sample with the lowest number of reads. Any user-specified value must be lower than the number of reads in the sample with the lowest number of reads.
}
\item{rarefy_each_step}{
Whether or not to rarefy each bootstrap sample separately. The default is TRUE (i.e., perform a new rarefaction on the OTU table for each bootstrap sample). However, for cases with large numbers of samples, large numbers of reads per sample or high microbial diversity, this can be slow; thus, there is the option to rarefy the microbiome samples a single time prior to bootstrapping.
}
}
\details{
\code{FourHbootstrap} generates bootstrapped samples of the 4H index. Each bootstrap begins by rarefying the OTU table such that all samples have the same library size. Once the OTU table has been rarefied, a subset, \code{sample_no}, of hosts is selected from each host class. Hosts can be sampled with or without replacement. The default is to sample hosts without replacement, such that all hosts in any given bootstrap sample are unique. Subsampled individuals from each host class are then used to determine core microbiomes. This allows calculation of the four dimensions of the 4H-index:
\deqn{\mathcal{U}=\frac{|\neg (P_1 \cap\ P_2) \cap\ H|}{|P_1 \cup\ P_2 \cup\ H}|}
\deqn{\mathcal{I}=\frac{|(P_1 \cap\ P_2) \cap\ H|}{|P_1 \cup\ P_2 \cup\ H}|}
\deqn{\mathcal{G}=\frac{|(P_1 \cup\ P_2) \cap\ H|}{|P_1 \cup\ P_2 \cup\ H}|}
\deqn{\mathcal{L} = 1 - \mathcal{U} -\mathcal{I} -\mathcal{G}}
where \eqn{P_1} is the set of core microbial taxa from the first parent species, \eqn{P_2} is the set of core microbial taxa from the second parent species, \eqn{H} is the set of core microbial taxa from hybrid organisms, and \eqn{|S|} denotes the cardinality of set \eqn{S}. \eqn{\mathcal{U}}, \eqn{\mathcal{I}}, \eqn{\mathcal{G}} and \eqn{\mathcal{L}} represent the extent of the Union, Intersection, Gain and Loss Models respectively (see Camper et al., 2023). A 4H-index is calculated for each bootstrap sample.

In addition to the four dimensions of the 4H-index, \code{FourHbootstrap} also calculates \eqn{\sigma} for each bootstrap sample:
\deqn{\sigma=\frac{P_1 \cap\ P_2}{P_1 \cup\ P_2}}
where \eqn{\sigma} is the fraction of parental microbial taxa that are found on both parent species. This will be used later for calculating the null plane of a given hybrid complex (see \code{FourHnullplane}).
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

FourHbootstrap(enterotype,hybrid_status,0.5,5,10)

}
