\name{averagetable}

\alias{averagetable}

\title{Computes a (products,descriptors) matrix}

\description{
  Returns the (products,descriptors) matrix with entries the means over panelists and sessions. \cr
  Computes analyses of variance automatically for a given model and a set of quantitative variables.
  Returns a data matrix where each row is associated with each
  category of a given categorical variable (in most cases, the categorical variable is the 
  \emph{product} variable), each column is associated with a quantitative variable, and each cell is
  the corresponding adjusted mean or mean.\cr
  Computes the average data table with respect to a categorical variable and a set 
  of quantitative variables.
 
            }

\usage{
averagetable(donnee, formul, subset = NULL, method = "coeff", 
    firstvar, lastvar = ncol(donnee), file = NULL)
}

\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{formul}{the model with respect to which the factor levels of the categorical variable of interest are calculated}  
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process}
  \item{method}{two possibilities, "coeff" (by default) or "mean"}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{file}{the name of the output file (by default, NULL and results are not in a file)}
}

\details{
  The \code{formul} parameter can be filled in for a given analysis of variance model.
  The \code{formul} parameter must begin with the categorical variable of interest (generally the \emph{product} variable)
  followed by the different other factors (and eventually their interactions) of interest. Classicially, one can used 
  \code{formul = "~Product+Panelist+Product:Panelist"}.
  In practise and in our type of applications, this function is very useful to obtain a data matrix 
  in which rows represent products and columns represent sensory descriptors. \cr
  If "mean" is assigned to the \code{method} parameter, then the \code{formul} parameter
  can be restricted to the sole variable of interest (generally the \emph{product} variable). \cr
  If data are balanced, the two options "mean" and "coeff" give the same results.
  
}

\value{
Return a matrix of dimension (\emph{p,q}), where \emph{p} is the number of categories of the qualitative variable 
of interest (in most cases, \emph{p} is the number of products)
and \emph{q} is the number of (sensory) descriptors. If "coeff" is assigned to the 
\code{method} parameter then the function \emph{averagetable} returns the matrix
of the adjusted means; if "mean" is assigned to the \code{method} parameter 
then the function averagetable returns the matrix of the means per category.
}

\references{P. Lea, T. Naes, M. Rodbotten. \emph{Analysis of variance for sensory data}. \cr
H. Sahai, M. I. Ageel. \emph{The analysis of variance}.
}

\author{François Husson \email{François.Husson@agrocampus-rennes.fr}}

\seealso{ \code{\link{aov}}}

\examples{
data(chocolates)
resaverage<-averagetable(sensochoc, formul = "~Product+Panelist",
    firstvar = 5)
coltable(magicsort(resaverage), level.upper = 6,level.lower = 4,
    main.title = "Average by chocolate")

res.pca = PCA(resaverage, scale.unit = TRUE)
}

\keyword{models}
