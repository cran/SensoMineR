\name{decat}

\alias{decat}

\title{DEscription of CATegories}

\description{
  This function is designed to point out the variables that are the most characteristic according 
  to the set of products in its whole, and to each of the products in particular. \cr

  This function is designed to test the main effect of a categorical variable (F-test) 
  and the significance of its coefficients (T-test) for a set of endogenous variables 
  and a given analysis of variance model.
  In most cases, the main effect is the product effect and the endogenous variables are
  the sensory descriptors.
}

\usage{
decat(donnee, formul, firstvar, lastvar = length(colnames(donnee)), 
    proba = 0.05, graph = TRUE)
}

\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{formul}{the model that is to be tested}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{proba}{the significance threshold considered for the analyses of variance (by default 0.05)}
  \item{graph}{a boolean, if TRUE a barplot of the P-values associated with the F-test of the product effet is displayed}
}

\details{
  The \code{formul} parameter must be filled in by an analysis of variance model and must 
  begin with the categorical variable of interest (e.g. the product effect)
  followed by the different other factors of interest (and their combinations).
  E.g.:\code{formul = "~Product+Panelist+Session"}.
}

\value{
A list containing the following elements:
  \item{tabF}{the V-test and the P-value of the F-test for each descriptor resulting from the analysis of variance model}
  \item{tabT}{a (products,descriptors) data frame, in which each cell is the Vtest for a given product and
  a given descriptor}
  \item{coeff}{a (products,descriptors) data frame, in which each cell is the coefficient resulting
  from the analysis of variance model for a given product and a given descriptor}
  \item{resF}{the V-test and the P-value for each descriptor resulting from the analysis of variance model, sorted in ascending order}
  \item{resT}{a list which elements are data frames, one data frame per product: 
  the coefficient, the P-value and the Vtest for each significant descriptor resulting from the analysis of variance model, sorted in descending order}
A barplot of the P-values associated with the F-test of the product effet.
}

\references{P. Lea, T. Naes, M. Rodbotten. \emph{Analysis of variance for sensory data}.\cr
H. Sahai, M. I. Ageel. \emph{The analysis of variance}.
}

\author{François Husson}

\seealso{ \code{\link{aov}}}

\examples{
data(chocolates)
## model (AOV): " descriptor = product + panelist "
resdecat<-decat(chocolates, formul="~Product+Panelist", firstvar = 5)
coltable(magicsort(resdecat$tabT), level.lower = -1.96, 
    level.upper = 1.96, main.title="Products' description")
barrow(resdecat$tabT)
barrow(t(resdecat$tabT), numr = 3, numc = 3)
barrow(resdecat$coeff, color = "orange") 
}

\keyword{models}
