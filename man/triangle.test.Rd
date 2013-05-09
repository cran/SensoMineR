\name{triangle.test}
\alias{triangle.test}


\title{Make a Triangle test for a set of products}

\description{
Make a Triangle test for a set of products.
}

\usage{
triangle.test (design, answer, preference = NULL)
}

\arguments{
  \item{design}{a data.frame corresponding to the design use to make the Triangle test (typically the ouput of the function \code{triangle.design}}
  \item{answer}{a vector of the answers of all the panelists; all the answer should be "X", "Y" or "Z"}
  \item{preference}{a vector of the preference of the panelists; all the answer should be "X", "Y" or "Z" (by default, there preference are not taken into account)}
}
\details{
    Triangle test: panelists receive three coded samples. They are told that two of the sample are the same and one is different.
    Panelists are asked to identify the odd sample. 
}
\value{
Returns a list of matrices. Each matrix give the reult for all the pair of products:
\item{nb.comp}{a matrix with the number of comparisons done for each pair of products;}
\item{nb.ident}{a matrix with the number of panelists who indicate the odd product for each pair of products;}
\item{p.value}{a matrix with the p-value of the Triangle tests for each pair of products;}
\item{nb.recognition}{estimation of the panelists who really perceived the difference between two product, for each pair of product;}
\item{maxML}{Maximum Likelihood of the estimation of the number of panelists who really perceive the difference between the products;}
\item{confusion}{estimation of the percentage of panelists who do not perceived the difference between two product, for each pair of product;}
\item{minimum}{minimum of panelists who should detect the odd product to can say that panelists perceive the difference between the products, for each pair of products;}
\item{preference}{number of times that product of row i is prefered that product in column j for the panelists who find the odd product.}
}


\author{François Husson}

\seealso{\code{\link{triangle.pair.test}}, \code{\link{triangle.design}}}

\examples{
design = triangle.design(nbprod = 4, nbpanelist = 6, bypanelist = 3)
answer = c("X","Y","Y","X","Z","X","Y","X","Z",
    "X","X","Z","X","Y","X","Z","X","Y")
triangle.test (design, answer)  
}
\keyword{models}
