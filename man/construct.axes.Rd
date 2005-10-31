\name{construct.axes}

\alias{construct.axes}

\title{Coordinates of individuals and illustrative individuals for PCA or MFA}

\description{
This function is especially designed to be used in a sensory data analysis context.
Returns the coordinates of the products when performing either PCA or MFA and 
the coordinates of the "partial" products when performing MFA.
Returns also the panelists' coordinates when projected as illustrative rows onto the products' space.
Produces graphs of products and descriptors from the output of PCA or MFA.
}

\usage{
construct.axes(ktableau,coord=c(1,2),scale.unit=TRUE,centerbypanelist=FALSE,scalebypanelist=FALSE)
}

\arguments{
  \item{ktableau}{a ktab made up of at least two blocks (or sets) of variables: the first
  block must be a two qualitative variables block with the \emph{panelist} and 
  the \emph{product} variables, the others are sensory descriptors blocks used to perform an MFA
  (to perform a PCA only one additional block to the first one is needed)}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{scale.unit}{boolean, if \code{TRUE} the descriptors are scaled to unit variance}
  \item{centerbypanelist}{center the data by panelist before the construction of the axes}
  \item{scalebypanelist}{scale the data by panelist before the construction of the axes}
}

\details{

The input data set is an object of class \code{ktab} (i.e. a list of data frames with the same row.names in common), 
for which the first element of the list (i.e. the first data frame) must correspond to a 2
qualitative variables block, where the first variable refers to the \emph{panelist} variable 
and the second to the \emph{product} variable. \cr

The ouptut of this function is a list with one element when performing PCA and two elements 
when performing MFA. The first element is the data frame of the coordinates of the products according to 
the whole panel (Panelist=0) and to the panelists. 
The second element is the data frame of the coordinates of the "partial products"
according to the whole panel (Panelist=0) and to the panelists. \cr

This function is necessary when calculating confidence ellipses for products.
}

\value{
A list containing the following elements:
  \item{moyen}{the coordinates of the products with respect to the panel and to each panelists}
  \item{partiel}{the coordinates of the \emph{partial} products with respect to the panel and to each panelists}
Returns also a correlation circle as well as a graph of individuals
}

\references{
Escofier, B. and Pagès, J. (1990) 
\emph{Analyses factorielles simples et multiples: objectifs, méthodes et interprétation} 
Dunod, Paris. 1--267. \cr

Escofier, B. and Pagès, J. (1994) 
Multiple factor analysis (AFMULT package). 
\emph{Computational Statistics and Data Analysis}, \bold{18}, 121--140.
}

\author{François Husson}

\seealso{ \code{\link{mfasenso}}}

\examples{
\dontrun{
data(chocolates)

## PCA (two blocks of variables)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,14),tabnames=c("JP","Gr1"))
axe <- construct.axes(ktab.donnee,scale.unit=TRUE)
 
## MFA (three blocks of variables)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,6,4,4),tabnames=c("JP","A","F","T"))
axe <- construct.axes( ktab.donnee,scale.unit=TRUE)
}
}
\keyword{multivariate}
