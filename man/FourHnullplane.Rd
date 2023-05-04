\name{FourHnullplane}
\alias{FourHnullplane}
\title{
Null Plane
}
\description{
Calculates and draws the null plane on a quaternary (four dimensional barycentric, Aitchison Simplex) plot.
}
\usage{
FourHnullplane(boots,col='red')
}
\arguments{
  \item{boots}{(Required) Output matrix from the \code{FourHbootstrap} function
}
\item{col}{The color of the plane in the quaternary plot. The default color is red.}

}
\details{
\code{FourHnullplane} uses the average value of \eqn{\sigma= (P_1 \cap\ P_2)/(P_1 \cup\ P_2)} across all bootstrap samples to draw a null plane on an existing quaternary plot. For any given value of \eqn{\theta = \mathcal{G}+\mathcal{L}}, the null plane is defined as follows:
\deqn{\mathcal{U}_{null} = (1-\sigma)(1-\theta)}
\deqn{\mathcal{I}_{null} = \sigma(1-\theta)}
The primary purpose of the null plane is visualization. In particular, the null plane reflects the ratio of \eqn{\mathcal(I)/(\mathcal{U}+\mathcal{I})} that would be expected if hybrids were equally likely to acquire any microbial taxa found on one or both parent species. 4H-indices that lie closer to the \eqn{\mathcal{I}} vertex than the null plane indicate that hybrid organisms are more likely to retain microbial taxa shared by both parent species than microbial taxa found on only one parent species (and by corollary, hybrid organisms are more likely to lose microbial taxa found on only one parent species than they are to lose microbial taxa found on both parent species). 4H-indices that lie closer to the \eqn{\mathcal{U}} vertex than the null plane indicate that hybrid organisms are more likely to retain microbial taxa found on only one parent species than microbial taxa found on both parent species (and by corollary, hybrid organisms are more likely to lose microbial taxa found on both parent species than they are to lose microbial taxa found on only one parent species).
}
\value{
This function outputs a plane on an existing quaternary plot.
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
boot_samples<-FourHbootstrap(enterotype,hybrid_status, 0.5,5,10)

#Plot the null plane on an existing quaternary plot
#An existing quaternary plot should already be open
FourHnullplane(boot_samples,col='red')

}
