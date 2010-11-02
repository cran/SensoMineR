\name{print.fasnt}

\alias{print.fasnt}

\title{Print Factorial Approach for Sorting Napping Task data (FASNT) results}

\description{
Print Factorial Approach for Sorting Napping Task data (FASNT) results.
}

\usage{
\method{print}{fasnt}(x, file = NULL, sep = ";", \dots)
}

\arguments{
  \item{x}{an object of class fasnt}
  \item{file}{A connection, or a character string naming the file to print to.  If NULL (the default), the results are not printed in a file}
  \item{sep}{character string to insert between the objects to print (if the argument file is not NULL}
  \item{\dots}{further arguments passed to or from other methods}
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fasnt}}}

\examples{
\dontrun{
data(smoothies)
res.fasnt <- fasnt(smoothies, first="catego",graph=F)
print.fasnt(res.fasnt, file="c:/fasnt.csv", sep = ";")
}
}
\keyword{print}
