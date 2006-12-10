\name{scalebypanelist}

\alias{scalebypanelist}

\title{Scale by panelist}

\description{
Returns a data frame with entries the means of the products over the sessions for 
the whole panel and for each panelist.
}

\usage{
scalebypanelist(matrice, center = TRUE, scale = FALSE, col.p, 
    col.j, firstvar, lastvar = ncol(matrice), method = "coeff")
}

\arguments{
  \item{matrice}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{center}{boolean, if T scores given by panelists are centered}  
  \item{scale}{boolean, if T scores given by panelists are scaled to unit varaince}  
  \item{col.p}{the position of one categorical variables of interest (the \emph{product} variable)}
  \item{col.j}{the position of one categorical variables of interest (the \emph{panelist} variable)}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{method}{the method to replace the missing values: "average" or "coeff" 
  (coefficients of the \emph{product} variable in the anova model)}
}

\value{
Returns a data frame of dimension (\emph{p*(1+j),k+2}), where \emph{p} is the number
of products, \emph{j} the number of panelists, and \emph{k} the number of sensory 
descriptors (the first two variables correspond to the \emph{panelist} and the \emph{product} variables). 
This data frame contains the means of the products over the sessions for 
the whole panel and for each panelist (data may be scaled to unit variance or not, this 
parameter is set by users).
}

\author{François Husson}

\examples{
data(chocolates)
res=scalebypanelist(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
res
}

\keyword{manip}
