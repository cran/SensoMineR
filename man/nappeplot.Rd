\name{nappeplot}
\alias{nappeplot}

\title{Plot panelists' tableclothe}
\description{Plot panelists' tableclothe.}
\usage{
nappeplot(donnee, numr = 2, numc = 2, color = "blue", lim = c(60,40))
}

\arguments{
  \item{donnee}{a data frame of dimension (\emph{p,2j}),
  where \emph{p} represents the number of products and \emph{j} the number 
  of panelists}
  \item{numr}{the number of tableclothe per row (by default 2)}
  \item{numc}{the number of tableclothe per column (by default 2)}
  \item{color}{the color used to display the products}
  \item{lim}{the size of the tableclothe}  
}
\details{
The data used here refer to a specific experiment, where panelists are 
asked to position products on a tableclothe of dimension \code{lim}, by default (60,40).
}
\value{
Returns as many graphs as there are panelists, each graph represents products 
positioned by a given panelist on a tablecloth
}
\references{
Pagès J. (2005). Collection and analysis of perceived product inter-distances using multiple factor analysis; 
application to the study of ten white wines from the Loire Valley. \emph{Food Quality and Preference}. 16 (7) pp. 642-649.
}
\author{François Husson}

\seealso{\code{\link{napping}}, \code{\link{pmfa}}, \code{\link{indscal}}}
\examples{
\dontrun{
data(napping)
nappeplot(napping.don)
}
}
\keyword{multivariate}
