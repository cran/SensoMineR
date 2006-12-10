\name{boxprod}

\alias{boxprod}

\title{Boxplot per category with respect to a categorical variable and a set of quantitative variables}

\description{
  Returns as many boxplots as there are categories for a given categorical variable of interest (in most cases, the
  \emph{product} variable). The boxplots are automatically generated for all the quantitative variables 
  (in our type of applications, variables are often sensory descriptors).
}

\usage{
boxprod(donnee, col.p, firstvar, lastvar = ncol(donnee), 
    numr = 2, numc = 2)
}

\arguments{
  \item{donnee}{a data frame}
  \item{col.p}{the position of the categorical variable of interest}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{numr}{the number of boxplots per row (by default 2)}
  \item{numc}{the number of boxplots per column (by default 2)}
}

\details{
Missing values are ignored when forming boxplots. 
}

\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.\cr
Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey, P. A. (1983) \emph{Graphical Methods for Data Analysis}. Wadsworth & Brooks/Cole. 
}

\author{ 
François Husson \email{François.Husson@agrocampus-rennes.fr} \cr
Sébastien Lê \email{Sebastien.Le@agrocampus-rennes.fr}
}

\seealso{ 
\code{\link{boxplot}} which does the computation,
\code{\link{bxp}} for the plotting and more examples; and \code{\link{stripchart}} 
for an alternative (with small data sets).
}

\examples{
data(chocolates)
boxprod(sensochoc, col.p = 4, firstvar = 5, numr = 2, numc = 2)
}

\keyword{univar}
