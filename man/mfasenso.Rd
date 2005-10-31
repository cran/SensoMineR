\name{mfasenso}

\alias{mfasenso}

\title{Make MFA (or PCA) with illustrative individuals}

\description{

}

\usage{
mfasenso(ktab, ktab.illu = NULL, scale.unit = TRUE, nbcoord = 2, 
    poids = NULL)
}

\arguments{
  \item{ktab}{ktableau with active individuals}
  \item{ktab.illu}{ktableau with illustrative individuals}  
  \item{scale.unit}{Boolean, TRUE if the variables are scaled to unit}  
  \item{nbcoord}{number of coordinates}  
  \item{poids}{vector of the weight of the variables (by default, all the variables hace the same weight)}  
}

\details{
Factor analysis is done. Two ktable are required, the first one concern the active individuals,
the second one concerns the illustrative individuals. The outputs are the coordinates
of the active and illustrative individuals (the mean points and the partial points if the method used is an MFA).
If in the ktables, there is only one group of variable, then a PCA is done (else an MFA).
}

\value{

}

\references{
Escofier B., Pagès J. (1988, 1990, 1993, 1998) Analyses factorielles simples et multiples. Objectifs méthodes et interprétation. Dunod. Paris.
Pagès, J. & Husson, F. (2005) Multiple Factor Analysis with confidence ellipses: a useful methodology to study the relations between sensory and instrumental data. Journal of chemometrics.
}

\author{François Husson}

\examples{ 
## PCA 
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],
    chocolates[,4], chocolates[,-(1:4)]), blocks = c(2,14),
    tabnames = c("JP","Desc"))

## MFA 
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],
    chocolates[,4], chocolates[,-(1:4)]), blocks = c(2,6,8),
    tabnames = c("JP","A-F","T-S"))
}

\keyword{dplot}
\keyword{internal}
