\name{panelperf}

\alias{panelperf}

\title{Panel's performance according to its capabilities to dicriminate between products}

\description{
Computes automatically P-values associated with the F-test as well as the residual term for a given analysis of variance model.
}

\usage{
panelperf(donnee, formul, subset = NULL, firstvar, 
    lastvar = ncol(donnee), random = TRUE)
}

\arguments{
  \item{donnee}{a data frame}
  \item{formul}{the model that is to be tested}
  \item{subset}{cf. function \code{\link{lm}} in the \pkg{stats} package}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee}}
  \item{random}{boolean, effect should be possible as fixed or random (default as random)}
}

\details{
  The \code{formul} parameter must be filled in by an analysis of variance model and must 
  begin with the categorical variable of interest (e.g. the product effect)
  followed by the different other factors of interest (and their combinations).
  E.g.:\code{formul = "~Product+Session"}.
}

\value{ 
A list containing the following components:
  \item{p.value}{a matrix of dimension (\emph{k,m}) of \emph{P-values} associated with the F-test for the 
  \emph{k} descriptors and the \emph{m} factors and their combinations considered in the analysis
  of variance model of interest}
  \item{variability}{a matrix of dimension (\emph{k,m}) where the entries correspond to the percentages of variability
  due to the effects introduced in the analysis of variance model of interest}
  \item{res}{a vector of dimension \emph{k} of residual terms for the analysis
  of variance model of interest}
  \item{r2}{a vector of dimension \emph{k} of r-squared for the analysis of variance model of interest}
}

\references{P. Lea, T. Naes, M. Rodbotten. \emph{Analysis of variance for sensory data}.\cr
H. Sahai, M. I. Ageel. \emph{The analysis of variance}.
}

\author{François Husson, Sébastien Lê}

\seealso{ \code{\link{paneliperf}}, \code{\link{aov}}}

\examples{
data(chocolates)
res=panelperf(sensochoc, firstvar = 5, formul = "~Product+Panelist+
    Session+Product:Panelist+Session:Product+Panelist:Session")
## Sort results by product p.values.
coltable(magicsort(res$p.value, sort.mat = res$p.value[,1], bycol = FALSE,
    method = "median"), main.title = "Panel performance (sorted by product P-value)")
}

\keyword{models}
