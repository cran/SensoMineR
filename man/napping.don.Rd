\name{napping.don}

\alias{napping.don}

\docType{data}

\title{An example of Napping data}

\description{
The data used here refer to 10 different French wines evaluated by 11 panelists.
They were asked to position the wines on a tableclothe of dimension (60,40).
}

\usage{data(napping)}

\format{
 A data frame of dimension (10,22): each row represents
 a French wine, each couple (Xi,Yi) represents the 
 coordinates of the wines positioned on a tableclothe for a given panelist.
}

\source{
Département de mathématiques appliquées, Agrocampus Rennes
}

\examples{
\dontrun{
data(napping)
nappeplot(napping.don)
}
}
\keyword{datasets}
