\name{FourHpreanalysis}
\alias{FourHpreanalysis}
\title{
Bootstraps of the percentage of gamma diversity that is core
}
\description{
Calculates the percentage of the total gamma diversity that is part of the core microbiota for bootstrapped samples drawn from an object of class phyloseq containing host-associated microbial community data from both hybrid hosts and hosts of each parent species. Prints out a warning when the percentage of the total gamma diversity that is core for hybrid hosts is less than half of the average percentage that is core for the two parent populations.
}
\usage{
FourHpreanalysis(x, class_grouping, core_fraction, boot_no, sample_no,
core_average_abundance = 0, core_max_abundance = 0,
core_all_average_abundance = 0, core_all_max_abundance = 0,
replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL)
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
\item{core_average_abundance}{
The minimum average abundance at which a microbial taxon must be found across all hosts within a class to be considered part of a host class's 'core' microbiome. The default is zero (i.e., it does not need to occur at any minimum average abundance)
}
\item{core_max_abundance}{
The minimum abundance at which a microbial taxon must be found on at least one host within a class to be considered part of a host class's 'core' microbiome. The default is zero (i.e., it does not need to occur at any minimum abundance on any single host)
}
\item{core_all_average_abundance}{
The minimum average abundance at which a microbial taxon must be found across all hosts (regardless of class) to be considered part of a host class's 'core' microbiome. The default is zero (i.e., it does not need to occur at any minimum average abundance)
}
\item{core_all_max_abundance}{
The minimum abundance at which a microbial taxon must be found on at least one host (regardless of class) to be considered part of a host class's 'core' microbiome. The default is zero (i.e., it does not need to occur at any minimum abundance on any single host)
}
\item{replace_hosts}{
Whether or not to replace hosts during bootstrapping. The default is FALSE (i.e., all hosts in a given bootstrap are unique).
}
\item{reads}{
The number of reads to rarefy to for each microbiome sample. The default is to use the number of reads in the sample with the lowest number of reads. Any user-specified value must be lower than the number of reads in the sample with the lowest number of reads.
}
\item{rarefy_each_step}{
Whether or not to rarefy each bootstrap sample separately. The default is TRUE (i.e., perform a new rarefaction on the OTU table for each bootstrap sample). However, for cases with large numbers of samples, large numbers of reads per sample or high microbial diversity, this can be slow; thus, there is the option to rarefy the microbiota samples a single time prior to bootstrapping.
}
\item{seed}{
Optional seed for random number generator. The default is NULL (i.e., seed picked at random at time of simulation).
}
}
\details{
\code{FourHpreanalysis} generates bootstrapped samples of the percentage of gamma diversity for hybrids and each parent species that is part of the core. Each bootstrap begins by rarefying the OTU table such that all samples have the same library size. Once the OTU table has been rarefied, a subset, \code{sample_no}, of hosts is selected from each host class. Hosts can be sampled with or without replacement. The default is to sample hosts without replacement, such that all hosts in any given bootstrap sample are unique. Subsampled individuals from each host class are then used to determine overall gamma diversity and diversity of core microbiotas for hybrid hosts and hosts of each parent species. This allows calculation of the percentage of gamma diversity that is part of the core. When this percentage among hybrids is less than half of what it is, on average, among parent species, a warning is printed.
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

FourHpreanalysis(enterotype,hybrid_status,0.5,5,10)

}
