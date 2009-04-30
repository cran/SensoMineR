\name{chocolates}

\alias{chocolates}

\docType{data}

\title{Chocolates data}

\description{
The data used here refer to six varieties of chocolates sold in France.\cr

- For the sensory description: each chocolate was evaluated twice by 29 panelists according to 14
sensory descriptors;\cr
- For the hedonic data: each chocolate was evaluated on a structured scale from 0 to 10, 
by 222 consumers, according to their liking (0) or disliking (10);\cr
- For the sensory panels description: each chocolate was evaluated by 7 panels according to 14
sensory descriptors.
}

\format{

There are three data frames:
- sensochoc: a data frame with 348 rows and 19 columns: 5 qualitative variables
 (Panelist, Session, Form, Rank, Product) and 14 sensory descriptors;\cr
- hedochoc: a data frame with 6 rows and 222 columns: each row corresponds to a chocolate
 and each column to the hedonic scores given by one of the 222 consumers participating in the study;\cr
- sensopanels: a data frame with 6 rows and 98 columns: each row corresponds to a chocolate
 and each column to the mean over the panelists of a given panel 
 according to a sensory descriptor.
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
