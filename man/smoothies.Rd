\name{smoothies}

\alias{smoothies}

\docType{data}

\title{Smoothies}

\description{
The data used here refer to 8 smoothies on which 24 consumers performed a sorted napping task.

}

\usage{data(smoothies)}

\format{
A data frame with 8 rows (the number of smoothies) and 72 columns (the number of consumers * 3).
For each consumer, we have the coordinates of the products on the tablecloth associated with napping on the one hand and the partitionning variable associated with categorization on the other hand.
The columns are grouped by consumer.
}

\source{
Département de mathématiques appliquées, AGROCAMPUS OUEST Centre de Rennes
}

\examples{
\dontrun{
data(smoothies)
## Example of FASNT
res.fasnt<-fasnt(smoothies,first="nappe")
}
}

\keyword{datasets}