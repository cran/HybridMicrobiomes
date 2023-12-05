\name{FourHbootstrap}
\alias{FourHbootstrap}
\title{
Bootstraps of the incidence-based 4H Index
}
\description{
Calculates incidence-based 4H Indices for bootstrapped samples drawn from an object of class phyloseq containing host-associated microbial community data from both hybrid hosts and hosts of each parent species.
}
\usage{
FourHbootstrap(x, class_grouping, core_fraction, boot_no, sample_no,
core_average_abundance = 0, core_max_abundance = 0,
core_all_average_abundance = 0, core_all_max_abundance = 0,
replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,dist='Jaccard')
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
Whether or not to rarefy each bootstrap sample separately. The default is TRUE (i.e., perform a new rarefaction on the OTU table for each bootstrap sample). However, for cases with large numbers of samples, large numbers of reads per sample or high microbial diversity, this can be slow; thus, there is the option to rarefy the microbiome samples a single time prior to bootstrapping.
}
\item{seed}{
Optional seed for random number generator. The default is NULL (i.e., seed picked at random at time of simulation).
}
\item{dist}{
How to count taxa present on one, two or three host classes when calculating the 4H Index. Use 'Jaccard' to count each taxon once, regardless of the number of host classes it occurs on (i.e., inspired by the Jaccard Index for beta diversity). Use 'Sorensen' to count each taxon proportional to the number of host classes it occurs on (i.e., inspired by the Sorensen Index for beta diversity). The default is 'Jaccard'.
}
}
\details{
\code{FourHbootstrap} generates bootstrapped samples of the incidence-based 4H Index. Each bootstrap begins by rarefying the OTU table such that all samples have the same library size. Once the OTU table has been rarefied, a subset, \code{sample_no}, of hosts is selected from each host class. Hosts can be sampled with or without replacement. The default is to sample hosts without replacement, such that all hosts in any given bootstrap sample are unique. Subsampled individuals from each host class are then used to determine core microbiomes. This is done by selecting all microbial taxa that occur on a threshold number of individuals of that host class (alternately, different abundance thresholds can also be used, see above). Richness of core microbiomes shared among various host classes are then used to calculate a 4H Index for each bootstrap sample.  Specifically, we define:

\deqn{a=|(P_1 \cap P_2) \cap H|}
\deqn{b=|(P_1 \cup P_2) \cap H|-a}
\deqn{b_1=|P_1 \cap H|-a}
\deqn{b_2=|P_2 \cap H|-a}
\deqn{c=|H|-a-b}
\deqn{d=|P_1 \cup P_2 \cup H|-a-b-c}
\deqn{d_1=|P_1|-|P_1 \cap P_2| - |P_1 \cap H|+a}
\deqn{d_2=|P_2|-|P_1 \cap P_2| - |P_2 \cap H|+a}
\deqn{d_12=d-d_1-d_2}

where \eqn{P_1} is the set of core microbial taxa from the first parent species, \eqn{P_2} is the set of core microbial taxa from the second parent species, \eqn{H} is the set of core microbial taxa from hybrid organisms, and \eqn{|S|} denotes the cardinality of set \eqn{S}. Thus, \eqn{a} is the number of microbial taxa shared by both progenitors and the hybrid, \eqn{b} is the number of microbial taxa shared by one progenitor (but not both) and the hybrid, \eqn{b_1} is the number of microbial taxa shared by the first progenitor and the hybrid, \eqn{b_2} is the number of microbial taxa shared by the second progenitor and the hybrid, \eqn{c} is the number of microbial taxa found only on the hybrid, \eqn{d} is the number of microbial taxa found only on one or both progenitors, \eqn{d_1} is the number of microbial taxa found only on the first progenitor, \eqn{d_2} is the number of microbial taxa found only on the second progenitor, and \eqn{d_12} is the number of microbial taxa found on both progenitors.

When dist='Jaccard', the 4H Index is calculated as:

\deqn{\mathcal{U}=\frac{b}{a+b+c+d}}
\deqn{\mathcal{U}_1=\frac{b_1}{a+b+c+d}}
\deqn{\mathcal{U}_2=\frac{b_2}{a+b+c+d}}
\deqn{\mathcal{I}=\frac{a}{a+b+c+d}}
\deqn{\mathcal{G}=\frac{c}{a+b+c+d}}
\deqn{\mathcal{L}=\frac{d}{a+b+c+d}=1 - \mathcal{U} -\mathcal{I} -\mathcal{G}}

\eqn{\mathcal{U}}, \eqn{\mathcal{I}}, \eqn{\mathcal{G}} and \eqn{\mathcal{L}} represent the extent of the Union, Intersection, Gain and Loss Models respectively (see Camper et al., 2023) and \eqn{\mathcal{U}_1} and \eqn{\mathcal{U}_2} partition the Union model into components associated with the first and second progenitor.

When dist = 'Sorensen', the 4H Index is calculated as:

\deqn{\mathcal{U}=\frac{2b}{3a+2b+c+d_1+d_2+2d_12}}
\deqn{\mathcal{U}_1=\frac{2b_1}{3a+2b+c+d_1+d_2+2d_12}}
\deqn{\mathcal{U}_2=\frac{2b_2}{3a+2b+c+d_1+d_2+2d_12}}
\deqn{\mathcal{I}=\frac{3a}{3a+2b+c+d_1+d_2+2d_12}}
\deqn{\mathcal{G}=\frac{c}{3a+2b+c+d_1+d_2+2d_12}}
\deqn{\mathcal{L}=\frac{d_1+d_2+2d_12}{3a+2b+c+d_1+d_2+2d_12}=1 - \mathcal{U} -\mathcal{I} -\mathcal{G}}



In addition to the four dimensions of the 4H Index, \code{FourHbootstrap} also calculates \eqn{\sigma} for each bootstrap sample:

Jaccard:
\deqn{\sigma=\frac{a+d_12}{a+d_12+d_1+d_2}}

Sorensen:
\deqn{\sigma=\frac{2(a+d_12)}{2(a+d_12)+d_1+d_2}}

where \eqn{\sigma} is the fraction of parental microbial taxa that are found on both parent species. This will be used later for calculating the null plane of a given hybrid complex (see \code{FourHnullplane}).

}
\value{
This function returns a matrix, where rows are individual bootstrap samples and columns are the four dimensions of the 4H Index, as well as \eqn{\sigma}, \eqn{\mathcal{U}_1}, and \eqn{\mathcal{U}_2}.
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
