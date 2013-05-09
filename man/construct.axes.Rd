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
construct.axes(matrice, coord = c(1,2), scale.unit = TRUE, group = NULL,
      name.group = NULL, centerbypanelist = FALSE, scalebypanelist = FALSE,
      method = "coeff")
}

\arguments{
  \item{matrice}{a data.frame made up of at least two qualitative variables (the \emph{panelist} and 
  the \emph{product} variables), the others are sensory descriptors used to perform an MFA or a PCA if group = NULL}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{scale.unit}{boolean, if \code{TRUE} the descriptors are scaled to unit variance}
  \item{group}{the number of variables in each group of variables when multiple factor analysis is performed (by default this parameter equals NULL and a PCA is performed)}  
  \item{name.group}{the names of the groups of variables when mfa is performed (if \code{group} differs from NULL)}  
  \item{centerbypanelist}{center the data by panelist before the construction of the axes}
  \item{scalebypanelist}{scale the data by panelist before the construction of the axes}
  \item{method}{the method to replace the missing values: "average" or "coeff" 
  (coefficients of the \emph{product} variable in the anova model)}
}

\details{

The input data set is an object of class \code{data.frame}, 
for which the two first columns are qualitative variables (the first variable refers to the \emph{panelist} variable 
and the second to the \emph{product} variable) and the others are quantitative. \cr

The ouptut of this function is a list with one element when performing PCA and two elements 
when performing MFA. The first element is the data frame of the coordinates of the products according to 
the whole panel (Panelist=0) and to the panelists. 
The second element is the data frame of the coordinates of the "partial products"
according to the whole panel (Panelist=0) and to the panelists. \cr

This function is necessary when calculating confidence ellipses for products.
}

\value{
A list containing the following elements:
  \item{eig}{a matrix with the component of the factor analysis (in row) and the eigenvalues, the inertia and the cumulative inertia for each component}
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

\seealso{ \code{\link[FactoMineR]{MFA}}}

\examples{
## Example1: PCA 
data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
 
## Example2: MFA (two groups of variables)
data(chocolates) 
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, group = c(6,8), 
    name.group = c("A-F","T-S"),scale.unit = TRUE)
}
\keyword{multivariate}
