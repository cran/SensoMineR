\name{magicsort}
\alias{magicsort}

\title{Returns a sorted data matrix}
\description{
Sort the rows and columns of a matrix in a "magic" order or by ascending (or descending) mean or median or geometrical mean.
}
\usage{
magicsort(matrice, sort.mat = matrice, method = "magic",
    byrow = TRUE, bycol = TRUE, ascending = TRUE)
}

\arguments{
  \item{matrice}{a data matrix to sort}
  \item{sort.mat}{sort the rows and columns according to the result of the PCA made on this matrix (by default the \code{matrice})}
  \item{method}{four types of calculations, magic (\code{"magic"}), (\code{"median"}), arithmetical (\code{"mean"}) or geometrical (\code{"geo"}) mean (by default magic)}
  \item{byrow}{boolean, if TRUE then data are sorted over the rows}
  \item{bycol}{boolean, if TRUE then data are sorted over the columns}
  \item{ascending}{boolean, if TRUE then data are sorted in ascending order}
}
\details{
Very useful function to compare results.
}

\author{François Husson, Sébastien Lê}

\examples{
## Example 1
data(chocolates)
resdecat<-decat(sensochoc, formul = "~Product", firstvar = 5, 
    graph = FALSE)
coltable(magicsort(resdecat$tabT), level.lower = -1.96, 
    level.upper = 1.96, main.title = "Products' description")

## Example 2
data(chocolates)
resperf<-paneliperf(sensochoc, 
    formul = "~Product+Panelist+Product:Panelist",
    formul.j = "~Product", col.j = 1, firstvar = 5, lastvar = 12,
    synthesis = FALSE, graph = FALSE)
res.sort=magicsort(resperf$prob.ind, method = "median")
coltable(res.sort, main.title = "P-values of the F-test by panelist")
}
\keyword{manip}
