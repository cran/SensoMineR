\name{print.fast}

\alias{print.fast}

\title{Print Factorial Approach for Sorting Task data (FAST) results}

\description{
Print Factorial Approach for Sorting Task data (FAST) results.
}

\usage{
print.fast(x, file = NULL, sep = ";", \dots)
}

\arguments{
  \item{x}{an object of class fast}
  \item{file}{A connection, or a character string naming the file to print to.  If NULL (the default), the results are not printed in a file}
  \item{sep}{character string to insert between the objects to print (if the argument file is not NULL}
  \item{\dots}{further arguments passed to or from other methods}
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fast}}}

\examples{
\dontrun{
data(perfume)
res.fast <- fast(perfume,graph=FALSE)
print.fast(res.fast, file="c:/essai.csv", sep = ";")
}
}
\keyword{print}
