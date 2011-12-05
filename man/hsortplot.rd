\name{hsortplot}
\alias{hsortplot}

\title{Plot consumers' hierarchical sorting}
\description{Plot consumers' hierarchical sorting}
\usage{
hsortplot(don, group, numr = 2, numc = 2)
}

\arguments{
  \item{don}{a data frame with n rows (products) and p columns (nested partitions for all consumers)}
  \item{group}{a list indicating the number of levels (nested partitions) for each consumer}
  \item{numr}{the number of hierarchical sorting per row (by default 2)}
  \item{numc}{the number of hierarchical sorting per column (by default 2)}  
}
\details{
The data used here refer to a specific experiment, where children were asked to provide hierarchical sorting (several nested partitions) from 16 cards.
}

\value{
Returns as many graphs as there are consumers, each graph represents hierarchical sorting provided by a consumer}


\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{\code{\link{fahst}}}
\examples{
\dontrun{
data(cards)
group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
hsortplot(cards,group.cards)
}
}
\keyword{multivariate}
