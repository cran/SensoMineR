\name{triangle.design}
\alias{triangle.design}


\title{Construct a design for triangle tests}

\description{
Construct a design to make triangle tests.
}

\usage{
triangle.design (nbprod , nbpanelist, bypanelist = nbprod*(nbprod-1)/2,
   labprod=1:nbprod, labpanelist=1:nbpanelist)
}

\arguments{
  \item{nbprod}{number of products to compare}
  \item{nbpanelist}{number of panelists who make the triangle test}
  \item{bypanelist}{number of expermient that each panelist can done (by default each panelist make all the comparisons between the products}
  \item{labprod}{name of the products (by default, the product are coded from 1 to the number of products}
  \item{labpanelist}{name of the panelists (by default, the panelists are coded from 1 to the number of panelists}
}
\details{
Triangle test: panelists receive three coded samples. They are told that two of the sample are the same and one is different.
Panelists are asked to identify the odd sample.}
\value{
Returns an data.frame of dimension (\emph{t,3}), where \emph{t} is the number of experiments.
In column 1, 2 and 3 the product to test are given.
The product in column 1 is by coded "X", in column 2 is coded by "Y" and in column 3 is coded by "Z".
Panelist should start by product "X", then "Y" and then by "Z".
}


\author{François Husson}

\seealso{\code{\link{triangle.test}}, \code{\link{triangle.pair.test}}}

\examples{
##Example 1
design1 = triangle.design (nbprod = 4, nbpanelist = 8)

##Example 2
design2 = triangle.design(nbprod = 4, nbpanelist = 6, bypanelist = 3,
  labprod=c("prod1","prod2","prod3","prod4"),
  labpanelist=c("John","Audrey","Peter","Martina","James","Lisa"))
  
}
\keyword{models}
