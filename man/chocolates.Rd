\name{chocolates}

\alias{chocolates}

\docType{data}

\title{Chocolates data}

\description{
The data used here refer to six varieties of chocolates sold in France.
Each chocolate was evaluated twice by 29 panelists according to 14
sensory descriptors.
}

\usage{data(chocolates)}

\format{
 A data frame with 348 rows and 19 columns: 5 qualitative variables
 (Panelist, Session, Form, Rank, Product) and 14 sensory descriptors.
 
}

\source{
Département de mathématiques appliquées, Agrocampus Rennes
}

\examples{
data(chocolates)
decat(chocolates,formul="~Product+Panelist",firstvar=5,graph=FALSE)
}

\keyword{datasets}
