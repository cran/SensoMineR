\name{napping}

\alias{napping}

\docType{data}

\title{Napping data}

\description{
The data used here refer to 10 different French wines evaluated by 11 panelists.\cr
They were asked to position the wines on a tableclothe of dimension (60,40).
They were asked to describe each wine using their own word list.
}

\usage{data(napping)}

\format{
There are two data frames:
- napping.don: A data frame of dimension (10,22): each row represents
 a French wine, each couple (Xi,Yi) represents the 
 coordinates of the wines positioned on a tableclothe for a given panelist;\cr
- napping.words: A data frame of dimension (10,14): each row represents
 a French wine, each column an attribute, each cell the number of times
 a given attribute was quoted for a given wine.
}

\source{
Département de mathématiques appliquées, Agrocampus Rennes
}

\examples{
\dontrun{
data(napping)
nappeplot(napping.don)
x11()
pmfa(napping.don, napping.words)
}
}
\keyword{datasets}
