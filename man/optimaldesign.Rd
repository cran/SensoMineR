\name{optimaldesign}
\alias{optimaldesign}

\title{ Construction of an optimal design }
\description{
Construction of an optimal design balanced for first order of carry-over effect.
}
\usage{
optimaldesign(nbPanelist, nbProd, nbProdByPanelist=nbProd, seed=NULL)
}

\arguments{
  \item{nbPanelist}{ Maximum number of panelists }
  \item{nbProd}{ Number of products }
  \item{nbProdByPanelist}{ Number of products that each panelist will evaluate}
  \item{seed}{ initialization of the algorithm }
}

\value{
List with
  \item{design }{design with the products evaluated by each panelist}
  \item{rank }{contingency table with the number of times each product is seen at each rank}
  \item{succ }{contingency table with the number of times each product follow each product}

}

\author{F. Husson}

\examples{
\dontrun{
optimaldesign(nbPanelist=10,nbProd=5,nbProdByPanelist=3)
}
}
\keyword{ math }
