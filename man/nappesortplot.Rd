\name{nappesortplot}
\alias{nappesortplot}

\title{Plot consumers' sorted tablecloth}
\description{Plot consumers' sorted tablecloth.}
\usage{
nappesortplot(donnee,first="nappe", numr = 2, numc = 2, lim = c(60,40))
}

\arguments{
  \item{donnee}{a data frame of dimension (\emph{p,3j}),
  where \emph{p} represents the number of products and \emph{j} the number 
  of consumers}
  \item{first}{2 possibilities: "nappe" if the napping variables first appear for each consumer or "catego" if it is the categorization variable}
  \item{numr}{the number of tablecloth per row (by default 2)}
  \item{numc}{the number of tablecloth per column (by default 2)}
  \item{lim}{the size of the tablecloth}  
}
\details{
The data used here refer to a specific experiment, where consumers are 
asked to position products on a tablecloth of dimension \code{lim}, by default (60,40) and to categorize them.
}
\value{
Returns as many graphs as there are consumers, each graph represents products 
positioned by a given consumer on a tablecloth and colored according to the categorization
}
\references{
Pag\`es, J., L\^e, S., Cadoret, M. (2010) \emph{The Sorted Napping: a new holistic approach in sensory evaluation}. Journal of Sensory Studies\cr
}
\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{\code{\link{fasnt}}}
\examples{
\dontrun{
data(smoothies)
nappesortplot(smoothies,first="nappe")
}
}
\keyword{multivariate}
