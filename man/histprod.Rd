\name{histprod}

\alias{histprod}

\title{Histogram for each descriptor}

\description{
Computes automatically histograms for a set of quantitative variables.
}

\usage{
histprod(donnee, firstvar, lastvar = ncol(donnee), numr = 2, 
    numc = 2, adjust = 1)
}

\arguments{
  \item{donnee}{a data frame}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{numr}{the number of histograms per row (by default 2)}
  \item{numc}{the number of histograms per column (by default 2)}
  \item{adjust}{the bandwidth used is actually 'adjust*bw'. This makes it easy to specify values like "half the default" bandwidth.}
}

\details{
Displays histograms with a common Y-axis as well as the local estimator of the density for each descriptor, hence the
\code{adjust} parameter to fill in. Displays also the normal distribution with mean and variance the respective values
estimated for each descriptor.
}

\author{
Sébastien Lê
}

\seealso{ 
\code{\link{density}}, \code{\link{hist}}
}

\examples{
data(chocolates)
histprod(sensochoc, firstvar = 5, lastvar = 10)
}

\keyword{univar}
