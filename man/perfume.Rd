\name{perfume}

\alias{perfume}

\docType{data}

\title{Perfume}

\description{
The data used here refer to 12 luxury perfumes categorized by 30 consumers.

}

\usage{data(perfume)}

\format{
A data frame with 12 rows (the number of perfumes) and 30 columns (the number of consumers): a cell corresponds either to the number of the group to which the product belongs for the consumer, or, in the case of "qualified" categorization, to the sequence of words associted with the group to which the product belongs for the consumer.
}

\source{
Applied Mathematics Department, AGROCAMPUS OUEST Centre de Rennes
}

\examples{
data(perfume)

## Example of FAST
res.fast <- fast(perfume)

}

\keyword{datasets}