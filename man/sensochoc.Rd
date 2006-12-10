\name{sensochoc}

\alias{sensochoc}

\docType{data}

\title{Sensory data for 6 chocolates}

\description{
The data used here refer to the sensory description of six varieties of chocolates sold in France: 
each chocolate was evaluated twice by 29 panelists according to 14
sensory descriptors.
}

\format{

A data frame with 348 rows and 19 columns: 5 qualitative variables
 (Panelist, Session, Form, Rank, Product) and 14 sensory descriptors.
}

\usage{data(chocolates)}

\source{
Département de mathématiques appliquées, Agrocampus Rennes
}

\examples{
data(chocolates)
decat(sensochoc, formul = "~Product+Panelist", firstvar = 5, graph = FALSE)
}

\keyword{datasets}
