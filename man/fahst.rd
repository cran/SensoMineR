\name{fahst}

\alias{fahst}

\title{Factorial Approach for Hierarchical Sorting Task data}

\description{
Perform Factorial Approach for Hierarchical Sorting Task data (FAHST) on a table where the rows (i) are products and the columns (j) are for each consumer the partitionning variables associated with nested sorting. The columns are grouped by consumer.
For the partitionning variables, the label associated with a group can be an arbirary label (for example G1 for group 1, \emph{etc.}) or the words associated with the group in the case of qualified hierarchical sorting.}


\usage{fahst(don,group,alpha=0.05,graph=TRUE,axes=c(1,2),name.group=NULL,ncp=5,B=200,ncp.boot=2)}

\arguments{
  \item{don}{a data frame with n rows (products) and p columns (nested partitions for all consumers)}
  \item{group}{a list indicating the number of levels (nested partitions) for each consumer}
  \item{alpha}{the confidence level of the ellipses}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{name.group}{a vector containing the name of the consumers (by default, NULL and the consumers are named J1, J2 and so on)}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{B}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{ncp.boot}{number of dimensions used for the Procrustean rotations to build confidence ellipses (by default 2)}
}

\value{
A list containing the following elements:
  \item{eig}{a matrix containing all the eigenvalues, the percentage of variance and the cumulative percentage of variance}
  \item{ind}{a list of matrices containing all the results for the products (coordinates, square cosine, contributions)}
  \item{var}{a list of matrices containing all the results for the categories of the different nested partitions (coordinates, square cosine, contributions, v.test)}
  \item{group}{a list of matrices containing all the results for consumers (coordinates, square cosine, contributions)}
  \item{call}{a list with some statistics}
}

\references{
Cadoret, M., Le, S., Pages, J. (2010) \emph{A new approach for analyzing hierarchical sorting task data}. Sensometrics conference. Rotterdam, the Netherlands\cr
}

\author{Marine Cadoret, Sebastien Le \email{sebastien.le@institut-agro.fr}}

\examples{
\dontrun{
data(cards)
## Example of FAHST results
group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
res.fahst<-fahst(cards,group=group.cards)
}
}
\keyword{multivariate}