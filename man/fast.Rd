\name{fast}

\alias{fast}

\title{Factorial Approach for Sorting Task data}

\description{
Perform Factorial Approach for Sorting Task data (FAST) on a table where the rows (i) are products and the columns (j) are consumers.
A cell (i,j) corresponds either to the number of the group to which the product i belongs to for the consumer j, or, in the case of "qualified" categorization, to the sequence of words associted to the group of which the product belongs to (i) for the consumer j.
}

\usage{fast(don,alpha=0.05,mot_min=2,graph=TRUE,ncp=5,B=200)}

\arguments{
  \item{don}{a data frame with n rows (products) and p columns (assesor : categorical variables)}  
  \item{alpha}{the confidence level of the ellipses}
  \item{mot_min}{minimum sample size for the word selection in textual analysis}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{B}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
}

\value{
A list containing the following elements:
  \item{eig}{a matrix containing all the eigenvalues, the percentage of variance and the cumulative percentage of variance}
  \item{var}{a list of matrices containing all the results for the categories (coordinates, square cosine, contributions, v.test)}
  \item{ind}{a list of matrices containing all the results for the products (coordinates, square cosine, contributions)}
  \item{group}{a list of matrices containing all the results for consumers (coordinates, square cosine, contributions)}
  \item{cooccur}{the reordered co-occurrence matrix among products}
  \item{reord}{the reordered matrix products*consumers}
  \item{cramer}{the Cramer's V matrix between all the consumers}
  \item{call}{a list with some statistics}
}

\references{
Cadoret, M., L\^e, S., Pag\`es, J. (2008) \emph{A novel Factorial Approach for analysing Sorting Task data}. 9th Sensometrics meeting. St Catharines, Canada\cr
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\examples{
\dontrun{
data(perfume)
## Example of fast results
res.fast<-fast(perfume)
}
}
\keyword{datasets}