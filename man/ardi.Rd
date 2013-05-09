\name{ardi}
\alias{ardi}


\title{Automatic Research of DIvergences between scores}

\description{
Spot the most singular or particular data with respect to all descriptors and to two
qualitative variables and all their possible categories combinations.\cr
Computes the highest differences between all the 
categories of the variables \emph{product}, \emph{panelist} and all their possible combinations,
with respect to a set of quantitative variables (the sensory descriptors).
}

\usage{
ardi(donnee, col.p, col.j, firstvar, lastvar = ncol(donnee), 
      nbval = 10, center = TRUE, scale = FALSE)
}


\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{col.p}{the position of the \emph{product} variable}
  \item{col.j}{the position of the \emph{panelist} variable}
  \item{firstvar}{the position of the first sensory descriptor}
  \item{lastvar}{the position of the last sensory descriptor (by default the last column of \code{donnee})}
  \item{nbval}{the number of highest divergences to be displayed}
  \item{center}{by default, data are mean centered by panelist}
  \item{scale}{by default, data are not scaled by panelist}
}

\details{
Step 1 For each quantitative variable, means by all the possible combinations (panelist,product) are computed. \cr
Step 2 Then, data are mean centered and scaled to unit variance by descriptor and the divergence
corresponds to the absolute value of the entries. \cr
Step 3 Means on divergences are computed by products or by panelists and then sorted.
}

\value{
A list containing the following elements:
  \item{tab}{a data frame (descriptors are mean centered per panelist and scaled to unit variance)}
  \item{panelist}{a data frame, by default the 10 highest divergences between panelists according to the sensory
  descriptors}
  \item{product}{a data frame, by default the 10 highest divergences between products according to the sensory
  descriptors}
  \item{combination}{a data frame, by default the 10 highest divergences between panelists and products according to the sensory
  descriptors}
}


\author{
François Husson, Sébastien Lê
}

\seealso{\code{\link{decat}}}
\examples{
\dontrun{
data(chocolates)
ardi(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
}
  }
\keyword{univar}
