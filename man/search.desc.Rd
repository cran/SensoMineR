\name{search.desc}

\alias{search.desc}

\title{Search for discriminating descriptors}

\description{
This function is designed to select the significant descriptors in a data frame
}

\usage{
search.desc(matrice, col.j, col.p, firstvar, 
      lastvar = ncol(matrice), level = 0.5)
}

\arguments{
  \item{matrice}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{col.j}{the position of the categorical variable which make the variability, panelist for sensory studies. The value of \code{col.j} can also be NULL if no categorical variables make the variability.}
  \item{col.p}{the position of the categorical variable of interest, product for sensory studies}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee}}
  \item{level}{the threshold (P-value) below which variables are considered as discriminating for
  the following analysis of variance model: \code{descriptor=col.p+col.j} }  
}

\value{
Returns a data frame with all the qualitative variables and only discriminating variables
}


\author{François Husson}

\examples{
data(chocolates)
## In this example, all the descriptos are discriminated
interesting.desc <- search.desc(sensochoc, col.j = 1, col.p = 4, 
    firstvar = 5, level = 0.5)
}

\keyword{models}
