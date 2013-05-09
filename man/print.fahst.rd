\name{print.fahst}

\alias{print.fahst}

\title{Print Factorial Approach for Hierarchical Sorting Task data (FAHST) results}

\description{
Print Factorial Approach for Hierarchical Sorting Task data (FAHST) results.
}

\usage{
\method{print}{fahst}(x, file = NULL, sep = ";", \dots)
}

\arguments{
  \item{x}{an object of class fahst}
  \item{file}{A connection, or a character string naming the file to print to.  If NULL (the default), the results are not printed in a file}
  \item{sep}{character string to insert between the objects to print (if the argument file is not NULL}
  \item{\dots}{further arguments passed to or from other methods}
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fahst}}}

\examples{
\dontrun{
data(cards)
group<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
res.fast <- fahst(cards,group,graph=F)
print.fahst(res.fahst, file="c:/fahst.csv", sep = ";")
}
}
\keyword{print}
