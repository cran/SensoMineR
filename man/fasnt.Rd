\name{fasnt}

\alias{fasnt}

\title{Factorial Approach for Sorting Napping Task data}

\description{
Perform Factorial Approach for Sorting Napping Task data (FASNT) on a table where the rows (i) are products and the columns (j) are for each consumer the coordinates of the products on the tablecloth associated with napping on the one hand and the partitionning variable associated with categorization on the other hand. The columns are grouped by consumer.
For the partitionning variable, the label associated with a group can be an arbirary label (for example G1 for group 1, \emph{etc.}) or the words associated with the group in the case of qualified sorted napping.}


\usage{fasnt(don,first="nappe",alpha=0.05,sep.word=";",word.min=5,graph=TRUE,axes=c(1,2),
    ncp=5,name.group=NULL,B=100,val=FALSE, B.val=200)}
        
\arguments{
  \item{don}{a data frame with n rows (products) and p columns (assesor : categorical variables)}
  \item{first}{2 possibilities: "nappe" if the napping variables first appear for each consumer or "catego" if it is the categorization variable} 
  \item{alpha}{the confidence level of the ellipses}   
  \item{sep.word}{the word separator character in the case of qualified sorted napping}
  \item{word.min}{minimum sample size for the word selection in textual analysis}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{name.group}{a vector containing the name of the consumers (by default, NULL and the group are named J1, J2 and so on)}
  \item{B}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}  
  \item{val}{boolean, if TRUE elements of validity are calculating (it is time consuming)}
  \item{B.val}{the number of simulations used to obtain the elements of validity}
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
  \item{validity}{the elements of validity calculated for the first eigenvalue and the ellipses}
  \item{call}{a list with some statistics}
}

\references{
Pag\`es, J., L\^e, S., Cadoret, M. (2010) \emph{The Sorted Napping: a new holistic approach in sensory evaluation}. Journal of Sensory Studies\cr
Cadoret, M., L\^e, S., Pag\`es, J. (2009) \emph{Combining the best of two worlds, the "sorted napping"}. SPISE. Ho Chi Minh City, Vietnam\cr
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\examples{
\dontrun{
data(smoothies)
## Example of FASNT results
res.fasnt<-fasnt(smoothies,first="nappe")
}
}
\keyword{multivariate}