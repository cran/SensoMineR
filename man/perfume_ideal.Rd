\name{perfume_ideal}

\alias{perfume_ideal}

\docType{data}

\title{Perfume Ideal Data}

\description{
The data used here refer to the sensory description of twelve perfumes (2 were duplicated).\cr

Each perfume was evaluated once by 103 Dutch consumers and described on 21 attributes according to the Ideal Profile Method.\cr
Both perceived and ideal intensities were asked. In addition, the overall liking is asked.\cr
}

\format{

A data frame made of 103*14=1442 rows and 2 qualitative variables (set and product),
21*2 attributes (perceived and ideal intensities) and overall liking.
}

\usage{data(perfume_ideal)}

\source{
OP&P Product Research, Utrecht, The Netherlands
}

\examples{
\dontrun{
data(perfume_ideal)
decat(perfume_ideal, formul = "~product+user", firstvar = 3, graph = FALSE)
}
}
\keyword{datasets}
