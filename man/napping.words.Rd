\name{napping.words}

\alias{napping.words}

\docType{data}

\title{An example of "illustrative" variables to enhance results from Napping data}

\description{
The data used here refer to 10 different French wines evaluated by 11 panelists.
They were asked to describe each wine using their own word list.
}

\usage{data(napping)}

\format{
 A data frame of dimension (10,14): each row represents
 a French wine, each column an attribute, each cell the number of times
 a given attribute was quoted for a given wine
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
