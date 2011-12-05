\name{perfume_fcp}

\alias{perfume_fcp}

\docType{data}

\title{Perfume data obtained by free choice profiling}

\description{
The data used here refer to 12 luxury perfumes described by 6 experts.

}

\usage{data(perfume_fcp)}

\format{
A data frame with 12 rows (the number of perfumes) and 47 columns.
}

\examples{
data(perfume_fcp)
res <- fcp(perfume_fcp, group = c(12,7,7,7,6,8))
}

\keyword{datasets}