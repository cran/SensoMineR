\name{fasnt}

\alias{fasnt}

\title{Factorial Approach for Sorting Napping Task data}

\description{
Perform Factorial Approach for Sorting Napping Task data (FASNT) on a table where the rows (i) are products and the columns (j) are for each consumer the coordinates of the products on the tablecloth associated with napping on the one hand and the partitionning variable associated with categorization on the other hand. The columns are grouped by consumer.
For the partitionning variable, the label associated with a group can be an arbirary label (for example G1 for group 1, \emph{etc.}) or the words associated with the group in the case of qualified sorted napping.}


\usage{fasnt(don,first="nappe",B=100,axes=c(1,2),alpha=0.05,ncp=5,
     graph=TRUE,name.group=NULL,sep.word=" ",word.min=5,ncp.boot=2)}

\arguments{
  \item{don}{a data frame with n rows (products) and p columns (assesor : categorical variables)}
  \item{first}{2 possibilities: "nappe" if the napping variables first appear for each consumer or "catego" if it is the categorization variable}
  \item{B}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{alpha}{the confidence level of the ellipses}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{name.group}{a vector containing the name of the consumers (by default, NULL and the group are named J1, J2 and so on)}
  \item{sep.word}{the word separator character in the case of qualified sorted napping}
  \item{word.min}{minimum sample size for the word selection in textual analysis}
  \item{ncp.boot}{number of dimensions used for the Procrustean rotations to build confidence ellipses (by default 2)}
}

\value{
A list containing the following elements:
  \item{eig}{a matrix containing all the eigenvalues, the percentage of variance and the cumulative percentage of variance}
  \item{ind}{a list of matrices containing all the results for the products (coordinates, square cosine, contributions)}
  \item{quali.var}{a list of matrices containing all the results for the categories of categorization (coordinates, square cosine, contributions, v.test)}
  \item{quanti.var}{a list of matrices containing all the results for the napping (coordinates, square cosine, contributions, v.test)}
  \item{group}{a list of matrices containing all the results for consumers (coordinates, square cosine, contributions)}
  \item{indicator}{a list of matrices containing different indicators for napping and categorization}
  \item{textual}{the results of the textual analysis for the products}
  \item{call}{a list with some statistics}
}

\references{
Pag\`es, J., Le, S., Cadoret, M. (2010) \emph{The Sorted Napping: a new holistic approach in sensory evaluation}. Journal of Sensory Studies\cr
Cadoret, M., Le, S., Pages, J. (2009) \emph{Combining the best of two worlds, the "sorted napping"}. SPISE. Ho Chi Minh City, Vietnam\cr
}

\author{Marine Cadoret, Sebastien Le \email{sebastien.le@institut-agro.fr}}

\examples{
\dontrun{
data(smoothies)
## Example of FASNT results
res.fasnt<-fasnt(smoothies,first="nappe",sep.word=";")
}
}
\keyword{multivariate}