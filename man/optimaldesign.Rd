\name{optimaldesign}
\alias{optimaldesign}

\title{ Construction of an optimal design }
\description{
Construction of an optimal design balanced for first order of carry-over effect.
}
\usage{
optimaldesign(nbPanelist, nbProd, nbProdByPanelist, nbPanelistMin = nbPanelist,  ordre = TRUE, weight = 0.5, graine = Sys.time(), nbDesignProd = 10, nbDesignOrdre = 50, matEssImp = NA )
}

\arguments{
  \item{nbPanelist}{ Maximum number of panelists }
  \item{nbProd}{ Number of products }
  \item{nbProdByPanelist}{ Number of products that each panelist will evaluate}
  \item{nbPanelistMin}{ Minimum number of panelists who will evaluate the products }
  \item{ordre}{ Boolean, if TRUE the order of presentation of the product to the panelist is given }
  \item{weight}{ Importance of the rank and of the carry-over effect. From 0 to 1, if 0 the design will only take into account the carry-over effect (and not the rank effect).}
  \item{graine}{ initialization of the algorithm }
  \item{nbDesignProd}{ Number of iteration of the algorithm to affect the products to the panelists }
  \item{nbDesignOrdre}{ Number of iteration of the algorithm for the rank of presentation}
  \item{matEssImp}{ Matrix of the imposed experiments}
}

\value{
List with
  \item{design }{Design with the products evaluated by each panelist}
  \item{rank }{Design with the products evaluated by each panelist and with the rank}

}

\author{E. Périnel, O. Tran, J. Mazet}
\references{
Périnel E. & Pagès J. (2003). Optimal nested cross-over designs in sensory analysis. \emph{Food Quality and Preference}.  15 (5). pp. 439-446.
}

\examples{
\dontrun{
optimaldesign(nbPanelist=10,nbPanelistMin=8,nbProd=5,nbProdByPanelist=3)
}
}
\keyword{ math }
