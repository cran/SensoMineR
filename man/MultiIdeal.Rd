\name{MultiIdeal}

\alias{MultiIdeal}

\title{Single vs. Multiple Ideal}

\description{
By the use of confidence ellipses, this procedure checks whether consumers
associate the different products tested to a single or to multiple ideals.
}

\usage{
MultiIdeal(dataset, col.p, col.j, id.recogn, level.search.desc=0.2, correct=FALSE,
            nbchoix=NULL, nbsimul=500, coord=c(1,2))
}

\arguments{
  \item{dataset}{A matrix with at least two qualitative variables (\emph{consumer} and \emph{products})
  and a set of quantitative variables containing at least 2*A variables (for both \emph{perceived} and \emph{ideal} intensities)}
  \item{col.p}{The position of the \emph{product} variable }
  \item{col.j}{The position of the \emph{consumer} variable }
  \item{id.recogn}{The sequence in the variable names which distinguish the ideal
  variables from the sensory variables. This sequence should be fixed and unique. \cr
  Each ideal variable should be preceeded by the corresponding perceived intensity variable. }
  \item{level.search.desc}{the threshold above which a descriptor is not considered as discriminant according to AOV model "descriptor=Product+Panelist".}
  \item{correct}{Boolean, define whether the ideal products should be corrected from the difference in the use of the scale or not}
  \item{nbchoix}{The number of consumers forming a virtual panel, by default the number of panelists in the original panel}
  \item{nbsimul}{The number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{coord}{A length 2 vector specifying the components to plot}
}

\details{
The procedure of MultiIdeal, step by step:\cr
Step 1: the sensory and ideal variables are separated into two tables.\cr
Step 2: the product space is created by PCA on the averaged sensory table (averaged by product).\cr
Step 3: the ideal information (Product x Consumer) is projected as supplementary entities in this space.\cr
Step 4: confidence ellipses are created around the averaged ideal points associated to each product (using the consumer variability).\cr
}

\value{
Returns a matrix with the P-values of the Hotelling's T2 tests for each pair of products.
}

\references{Worch, T., & Ennis, J.M. (2013).
Investigating the single ideal assumption using Ideal Profile Method. 
\emph{Food Quality and Preference}.}

\author{Worch Thierry (thierry@qistatistics.co.uk)}

\seealso{ \code{\link{panellipse}}}

\examples{
\dontrun{
data(perfume_ideal)
res <- MultiIdeal(perfume_ideal, col.p=2, col.j=1, id.recogn="id_", 
    level.search.desc=0.2, nbsimul=500, coord=c(1,2))

# To run the analysis with all the attributes
res <- MultiIdeal(perfume_ideal, col.p=2, col.j=1, id.recogn="id_", 
    level.search.desc=1, nbsimul=500, coord=c(1,2))
}
}

\keyword{models}
