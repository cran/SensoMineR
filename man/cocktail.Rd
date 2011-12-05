\name{cocktail}

\alias{cocktail}

\docType{data}

\title{Cocktail data}

\description{
The data used here refer to 16 cocktails.

There are 3 files corresponding to the composition of the cocktails; the sensory description of the cocktails; the hedonic scores.

- For the composition of the cocktails: The mango, banana, orange and lemon concentration are known;\cr
- For the sensory description: each cocktail was evaluated by 12 panelists according to 13
sensory descriptors (only the average of each cocktail are given).
- For the hedonic data: each cocktail was evaluated on a structured scale from 0 to 10, 
by 100 consumers, according to their liking (0) or disliking (10).
}

\usage{data(cocktail)}

\format{
There are three data frames:
- compo.cocktail: a data frame with 16 rows and 4 columns: the composition of each cocktail is given for the 4 ingredients;\cr
- senso.cocktail: a data frame with 16 rows and 13 columns: each cocktail was evaluated by 12 panelists according to 13
sensory descriptors;\cr
hedo.cocktail: a data frame with 16 rows and 100 columns: each cocktail was evaluated on a structured scale from 0 to 10, 
by 100 consumers, according to their liking (0) or disliking (10).
}

\source{
Applied Mathematics Department, Agrocampus Rennes
}

\examples{
data(cocktail)
}

\keyword{datasets}
