\name{triangle.pair.test}
\alias{triangle.pair.test}


\title{Make a Triangle test for two products}

\description{
Make a Triangle test for two products.
}

\usage{
triangle.pair.test (nb.good, nb.answer)
}

\arguments{
  \item{nb.good}{number of panelists who identify the odd sample}
  \item{nb.answer}{number of panelists who make the triangle test}
}
\details{
    Triangle test: panelists receive three coded samples. They are told that two of the sample are the same and one is different.
    Panelists are asked to identify the odd sample. 
}
\value{
Returns 
\item{p.value}{the p-value of the Triangle test;}
\item{Estimation}{estimation by Maximum Likelihood of the number of panelists who really perceive the difference between the products;}
\item{ML}{Maximum Likelihood of the estimation of the number of panelists who really perceive the difference between the products;}
\item{minimum}{minimum of panelists who should detect the odd product to can say that panelists perceive the difference between the products.}
}


\author{François Husson}

\seealso{\code{\link{triangle.test}}, \code{\link{triangle.design}}}

\examples{
triangle.pair.test (11, 20)  
}
\keyword{models}
