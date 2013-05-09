\name{cards}

\alias{cards}

\docType{data}

\title{Cards}

\description{
The data used here refer to 16 cards (images) on which 30 children performed a hierarchical sorting task.

}

\usage{data(cards)}

\format{
A data frame with 16 rows (the number of cards) and 81 columns (the total number of levels provided by all children).
For each child, we have several qualitative variables corresponding to nested partitions: a partition corresponds to a level provided by the child.
The columns are grouped by child.
}

\source{
Département de mathématiques appliquées, AGROCAMPUS OUEST Centre de Rennes
}

\examples{
\dontrun{
data(cards)
## Example of FAHST
group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
res.fahst<-fahst(cards,group=group.cards)
}
}

\keyword{datasets}