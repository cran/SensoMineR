\name{cream_signa}

\alias{cream_signa}

\docType{data}

\title{Data description of the consumers who made the Ideal for the cream}

\description{
The data used here refer to the sensory description of 9 dessert chocolate creams.\cr

Each cream was evaluated once by 86 French consumers and described on 13 attributes according to the Ideal Profile Method.\cr
Both perceived and ideal intensities were asked. In addition, the overall liking is asked.\cr
}

\format{

A data frame made of 86*9=774 rows and 2 qualitative variables (panelist and product),
13*2 attributes (perceived and ideal intensities) and overall liking.
}

\usage{data(cream_signa)}

\source{Institut Agro, Melodie Sanchez, Sarah Sanchez
}

\examples{
\dontrun{
data(cream_signa)
data(cream_signa)
res.idmap <- IdMapConsumer(cream_id, cream_signa, col.p=2, col.j=1, col.lik=29, 
num.col.var.signa=c(1:12), conf.level=0.90, id.recogn="id_", color = FALSE, simusigni = 500) 

}
}
\keyword{datasets}
