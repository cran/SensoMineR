\name{afmult}
\alias{afmult}
\title{Multiple Factor Analysis}
\description{
Performs a multiple factor analysis, using an object of class \code{ktab}.
}
\usage{
afmult(X, option = c("lambda1", "inertia", "uniform", "internal"),
     coord = c(1,2), scannf = TRUE, nf = 3, cex = 0.8, 
     col = "steelblue4", font = 2, clabel = 0.8, scale.unit = TRUE)
}
\arguments{
  \item{X}{an object of class \code{ktab}}
  \item{option}{a string of characters for the weighting of arrays options:
    \item{\code{lambda1}}{weighting of group k by the inverse of the first eigenvalue of the k analysis}
    \item{\code{inertia}}{weighting of group k by the inverse of the total inertia of the array k}
    \item{\code{uniform}}{uniform weighting of groups}
    \item{\code{internal}}{weighting included in \code{X$tabw}}
  }
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{scannf}{a logical value indicating whether the eigenvalues bar plot should be displayed}
  \item{nf}{if scannf FALSE, an integer indicating the number of kept axes}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{col}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{font}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{clabel}{cf. the \pkg{ade4} package}
  \item{scale.unit}{a boolean, if TRUE (value set by default) then data are scaled to unit variance}
  }
\value{
  Returns a list including:
  \item{tab}{a data frame with the modified array}
  \item{rank}{a vector of ranks for the analyses}
  \item{eig}{a numeric vector with the all eigenvalues}
  \item{li}{a data frame with the coordinates of rows}
  \item{TL}{a data frame with the factors associated to the rows (indicators of table)}
  \item{co}{a data frame with the coordinates of columns}
  \item{TC}{a data frame with the factors associated to the columns (indicators of table)}
  \item{blo}{a vector indicating the number of variables for each table}
  \item{lisup}{a data frame with the projections of normalized scores of rows for each table}
  \item{link}{a data frame containing the projected inertia and the links between the arrays and the reference array}
}
\references{Escofier, B. and Pagès, J. (1994) Multiple factor analysis (AFMULT package), \emph{Computational Statistics and Data Analysis}, \bold{18}, 121--140.
}
\author{Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}}

\details{
The only difference between this function and the mfa programed by Daniel Chessel
lies in the calculations of the coordinates of the partial individuals.
}

\examples{
data(sensopanels)
resafmult<-afmult(ktab.data.frame(sensopanels, blocks = rep(14,7)),
    scannf = FALSE)
}

\keyword{multivariate}
