\name{pca}

\alias{pca}

\title{Principal components analysis}

\description{Performs a PCA and returns the individuals factor map and the variables factor map.} 

\usage{
pca(df, supind = NULL, supvar = NULL,
   row.w = rep(1,nrow(df)-length(supind))/(nrow(df)-length(supind)),
   scale.unit = TRUE, coord = c(1,2), graph = TRUE, 
   main.title = NULL, clabel = 1, cex = 0.7, font = 1, csub = 1, 
   col = "black", lty = 1)
}

\arguments{
  \item{df}{a data frame with \emph{n} rows (individuals) and \emph{p} columns (numeric variables)}
  \item{supind}{indexes of the illustrative individuals}
  \item{supvar}{indexes of the illustrative variables}
  \item{row.w}{an optional row weights (by default, uniform row weights)}
  \item{scale.unit}{a boolean, if TRUE (value set by default) then data are scaled to unit variance}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{main.title}{the title of the graph}
  \item{clabel}{if not NULL, a character size for the labels, used with  \code{par("cex")*clabel}}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{font}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{csub}{a character size for the legend, used with \code{par("cex")*csub}}
  \item{col}{color of the variables}
  \item{lty}{line type of the arrows}
}

\value{
Returns a list including:
  \item{eig}{a numeric vector with the all eigenvalues}
  \item{li}{a matrix with the coordinates of rows}
  \item{co}{a matrix with the coordinates of columns}
  \item{lisup}{a matrix with the coordinates of illustrative rows}
  \item{cosup}{a matrix with the coordinates of illustrative columns}

Returns the individuals factor map and the variables factor map.
}
  
}


\author{
François Husson, Sébastien Lê
}

\examples{
data(chocolates)
resaverage<-averagetable(chocolates, formul = "~Product+Panelist", 
    firstvar = 5)
pca(resaverage, scale.unit = TRUE)
}

\keyword{multivariate}
