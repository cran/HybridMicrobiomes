\name{FourHbootstrapA}
\alias{FourHbootstrapA}
\title{
Bootstraps of the abundance-based 4H Index
}
\description{
Calculates abundance-based 4H Indices for bootstrapped samples drawn from an object of class phyloseq containing host-associated microbial community data from both hybrid hosts and hosts of each parent species. The phyloseq object must, at a minimum, contain an OTU table AND sample data.
}
\usage{
FourHbootstrapA(x, class_grouping, core_fraction, boot_no, sample_no,
core_average_abundance = 0, core_max_abundance = 0,
core_all_average_abundance = 0, core_all_max_abundance = 0,
replace_hosts=FALSE, reads=NULL, rarefy_each_step=TRUE,seed=NULL,
dist='Bray-Curtis',representative='mean',rescale_core=FALSE, use_microViz = 'yes')
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
How to count taxa present on one, two or three host classes when calculating the 4H Index. Use 'Ruzicka' to count each taxon once, regardless of the number of host classes it occurs on (i.e., inspired by the Ruzicka Index for beta diversity). Use 'Bray-Curtis' to count each taxon proportional to the number of host classes it occurs on (i.e., inspired by the Bray-Curtis Index for beta diversity).
}
\item{representative}{
The method used to find representative microbial relative abundances for each host class. Options are 'mean', which is the average abundance of all microbial taxa found on individuals within the host class, and 'median' which is the median abundance of all microbial taxa found on individuals within the host class. If 'median' is chosen, then representative relative abundances are renormalized for each host class such that total abundances of all microbial taxa (core and non-core) are the same across host classes. The default is 'mean'.
}
\item{rescale_core}{
Whether or not to renormalize relative abundances across the subset of microbial taxa that comprise the core of each host class. Renormalizing the core microbial abundances of each host class results in the 4H Index being density invariant. However, this also constrains the 4H Index to a 2D plane, since there is no variation in total abundance across host classes. The default is to use raw relative abundances without renormalizing. Although this means that the 4H Index is not density invariant, in this case, differences in total abundance of the core microbiome between host classes are meaningful since they represent variation in the extent to which the overall microbiome of a host class is shared among individual animalss.
}
\item{use_microViz}{
Whether or not to use the microViz package to check your phyloseq object. If installation of the microViz package is problematic, this test can be skipped, but the user should verify themselves that a phyloseq object exists AND includes sample data (the sample data placeholder is necessary for the FourHbootstrapA code to work).
}

}
\details{
\code{FourHbootstrapA} generates bootstrapped samples of the abundance-based 4H Index. Each bootstrap begins by rarefying the OTU table such that all samples have the same library size. Once the OTU table has been rarefied, a subset, \code{sample_no}, of hosts is selected from each host class. Hosts can be sampled with or without replacement. The default is to sample hosts without replacement, such that all hosts in any given bootstrap sample are unique. Subsampled individuals from each host class are then used to (i) find representative relative abundances of each microbial taxon found on that host class and (ii) identify the core microbiome of the host class. Representative relative abundances are found by averaging or finding the median relative abundance of each microbial taxon across a particular host class. Core microbial taxa are identified by selecting all microbial taxa that occur on a threshold number of individuals of that host class (alternately, different abundance thresholds can also be used, see above). Once representative relative abundances of core microbial taxa have been determined, the core microbiome may or may not be renormalized (the default is to use the core microbial abundances without renormalization). Representative relative abundances are then used to determine whether certain host classes have an excess/deficit of abundance of any particular microbial taxa. Specifically, we define:

\deqn{A=\sum_{i=1}^{S}min(min(x_{iP_1},x_{iP_2}),x_{iH})=\sum_{i=1}^{S}\alpha_i}
\deqn{B=\sum_{i=1}^{S}[min(x_{iP_1}-\alpha_i,x_{iH}-\alpha_i)+min(x_{iP_2}-\alpha_i,x_{iH}-\alpha_i)]=\sum_{i=1}^{S}(\beta_{1i}+\beta_{2i})}
\deqn{B_1=\sum_{i=1}^{S}min(x_{iP_1}-\alpha_i,x_{iH}-\alpha_i)=\sum_{i=1}^{S}(\beta_{1i})}
\deqn{B_2=\sum_{i=1}^{S}min(x_{iP_2}-\alpha_i,x_{iH}-\alpha_i)=\sum_{i=1}^{S}(\beta_{2i})}
\deqn{C=\sum_{i=1}^{S}(x_{iH}-\alpha_i-\beta_{1i}-\beta_{2i})}
\deqn{D=\sum_{i=1}^{S}[x_{iP_1}+x_{iP_2}-\alpha_i-\beta_{1i}-\beta_{2i}-min(x_{iP_1},x_{iP_2})]}
\deqn{D_1=\sum_{i=1}^{S}[x_{iP_1}-\beta_{1i}-min(x_{iP_1},x_{iP_2})]}
\deqn{D_2=\sum_{i=1}^{S}[x_{iP_2}-\beta_{2i}-min(x_{iP_1},x_{iP_2})]}
\deqn{D_{12}=\sum_{i=1}^{S}[min(x_{iP_1},x_{iP_2})-\alpha_i]}

where \eqn{x_{iP_1}} is the relative abundance of microbial taxon \eqn{i} in a representative microbiome of the first parent species, \eqn{x_{iP_2}} is the relative abundance of microbial taxon \eqn{i} in a representative microbiome of the second parent species, \eqn{x_{iH}} is the relative abundance of microbial taxon \eqn{i} in a representative microbiome of the hybrid organism. Thus, \eqn{A} is the relative abundance of microbial taxa shared by both progenitors and the hybrid, \eqn{B} is the relative abundance of microbial taxa shared by one progenitor (but not both) and the hybrid, \eqn{B_1} is the relative abundance of microbial taxa shared by the first progenitor and the hybrid, \eqn{B_2} is the relative abundance of microbial taxa shared by the second progenitor and the hybrid, \eqn{C} is the relative abundance of microbial taxa found only on the hybrid, \eqn{D} is the relative abundance of microbial taxa found only on one or both progenitors, \eqn{D_1} is the relative abundance of microbial taxa found only on the first progenitor, \eqn{D_2} is the relative abundance of microbial taxa found only on the second progenitor, and \eqn{D_12} is the relative abundance of microbial taxa found on both progenitors.

When dist='Ruzicka', the 4H Index is calculated as:

\deqn{\mathcal{U}=\frac{B}{A+B+C+D}}
\deqn{\mathcal{U}_1=\frac{B_1}{A+B+C+D}}
\deqn{\mathcal{U}_2=\frac{B_2}{A+B+C+D}}
\deqn{\mathcal{I}=\frac{A}{A+B+C+D}}
\deqn{\mathcal{G}=\frac{C}{A+B+C+D}}
\deqn{\mathcal{L}=\frac{D}{A+B+C+D}=1 - \mathcal{U} -\mathcal{I} -\mathcal{G}}

\eqn{\mathcal{U}}, \eqn{\mathcal{I}}, \eqn{\mathcal{G}} and \eqn{\mathcal{L}} represent the extent of the Union, Intersection, Gain and Loss Models respectively (see Camper et al., 2023) and \eqn{\mathcal{U}_1} and \eqn{\mathcal{U}_2} partition the Union model into components associated with the first and second progenitor.

When dist = 'Bray-Curtis', the 4H Index is calculated as:

This allows calculation of the four dimensions of the 4H Index:
\deqn{\mathcal{U}=\frac{2B}{3A+2B+C+D_1+D_2+2D_{12}}}
\deqn{\mathcal{U}_1=\frac{2B_1}{3A+2B+C+D_1+D_2+2D_{12}}}
\deqn{\mathcal{U}_2=\frac{2B_2}{3A+2B+C+D_1+D_2+2D_{12}}}
\deqn{\mathcal{I}=\frac{3A}{3A+2B+C+D_1+D_2+2D_{12}}}
\deqn{\mathcal{G}=\frac{C}{3A+2B+C+D_1+D_2+2D_{12}}}
\deqn{\mathcal{L}=\frac{D_1+D_2+2D_{12}}{3A+2B+C+D_1+D_2+2D_{12}}=1 - \mathcal{U} -\mathcal{I} -\mathcal{G}}



In addition to the four dimensions of the 4H Index, \code{FourHbootstrapA} also calculates \eqn{\sigma} for each bootstrap sample:

Ruzicka:
\deqn{\sigma=\frac{A+D_{12}}{A+D_{12}+D_1+D_2}}

Bray-Curtis:
\deqn{\sigma=\frac{2(A+D_{12})}{2(A+D_{12})+D_1+D_2}}

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

FourHbootstrapA(enterotype,hybrid_status,0.5,5,10, use_microViz='no')

}
