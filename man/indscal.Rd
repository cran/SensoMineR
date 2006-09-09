\name{indscal}

\alias{indscal}

\title{Construct the Indscal model for Napping data type}

\description{
This version of the Indscal model is specially adapted to Napping data type, i.e. products (stimuli)
are positioned on a tableclothe by panelists, then their coordinates are used as input for the Indscal model.
}

\usage{
indscal(matrice, matrice.illu = NULL, maxit = 200, coord = c(1,2), 
    eps = 1/10^5)
}

\arguments{
  \item{matrice}{a data frame of dimension (\emph{p,2j}),
  where \emph{p} represents the number of products and \emph{j} the number 
  of panelists (two coordinates per panelist)}
  \item{matrice.illu}{a data frame with illustrative variables (with the same row.names in common as in \code{matrice})}  
  \item{maxit}{the maximum number of iterations until the algorithm stops}  
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{eps}{a threshold with respect to which the algorithm stops, i.e. when the difference between 
   the criterion function at step \emph{n} and \emph{n+1} is less than \code{eps}}  
}

\value{
Returns a list including:
  \item{W}{a matrix with the subject coordinates}
  \item{points}{a matrix with the stimuli (individuals) coordinates}
  \item{subvar}{a vector with the strain between each configuration and the stimuli configuration}
  \item{r2}{the strain criterion}

The functions returns the three following graphs: \cr
A stimuli representation, ie. a representation of the products\cr
A representation of the weights computed by the Indscal model.\cr
A correlation circle of the variables enhanced by illustrative variables (supplementary columns)
}

\references{
Carroll, J.D. & J.J. Chang (1970). Analysis of individual differences in multidimensional scaling via an N-way generalization of "Eckart-Young" decomposition. \emph{Psychometrika}, 35, 283-319.
}

\author{Peter Ellis \cr
François Husson}

\seealso{\code{\link{nappeplot}}, \code{\link{pmfa}}}

\examples{
\dontrun{
data(napping)
nappeplot(napping.don)
resindscal<- indscal(napping.don, napping.words)
x11()
prefpls(cbind(resindscal$points, napping.words))
x11()
pmfa(napping.don, napping.words, mean.conf = resindscal$points)
}
}
\keyword{multivariate}
