\name{interact}
\alias{interact}


\title{Estimation of interaction coefficients}

\description{
Computes automatically the interaction coefficients between two quantitative variables
\code{col.p} and \code{col.j} for the following model:
\code{"~col.p+col.j+col.p:col.j"}.
}

\usage{
interact(donnee, col.p, col.j, firstvar, lastvar = ncol(donnee))
}

\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{col.p}{the position of the \emph{product} effect for instance}
  \item{col.j}{the position of the \emph{panelist} effect for instance}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
}
\details{
In most cases \code{col.p} represents the \emph{product} effect, \code{col.j} represents the \emph{panelist} effect,
and the variables of interest are the sensory descriptors. The model considered is the following one:
\code{"~Product+Panelist+Product:Panelist"}.\cr
Data must be complete (but not necessarily balanced).
}
\value{
Returns an array of dimension (\emph{p,j,k}), where \emph{p} is the number of products, \emph{j} the number of panelists
and \emph{k} the number of sensory descriptors.
The entries of this array are the interaction coefficients between a panelist and a product for a given descriptor. \cr
For each sensory descriptor, returns a graph where each (panelist,product) interaction coefficient is displayed,
a graph where the contribution to the (panelist,product) interaction coefficient by product is displayed,
a graph where the contribution to the (panelist,product) interaction coefficient by panelist is displayed.
}


\author{François Husson}

\seealso{\code{\link{aov}}}

\examples{
\dontrun{
data(chocolates)
resinteract=interact(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
}
}
\keyword{models}
