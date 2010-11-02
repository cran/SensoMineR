\name{fast}

\alias{fast}

\title{Factorial Approach for Sorting Task data}

\description{
Perform Factorial Approach for Sorting Task data (FAST) on a table where the rows (i) are products and the columns (j) are consumers.
A cell (i,j) corresponds either to the number of the group to which the product i belongs for the consumer j, or, in the case of "qualified" categorization, to the sequence of words associted with the group to which the product i belongs for the consumer j.
}

\usage{fast(don,alpha=0.05,sep.words=";",word.min=5,graph=TRUE,axes=c(1,2),ncp=5,B=200,val=FALSE, B.val=200,label.miss=NULL)}

\arguments{
  \item{don}{a data frame with n rows (products) and p columns (assesor : categorical variables)}  
  \item{alpha}{the confidence level of the ellipses}
  \item{sep.words}{the word separator character in the case of qualified categorization}
  \item{word.min}{minimum sample size for the word selection in textual analysis}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{B}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{val}{boolean, if TRUE elements of validity are calculating (it is time consuming)}
  \item{B.val}{the number of simulations used to obtain the elements of validity}
  \item{label.miss}{label associated with missing groups in the case of incomplete data set}
}

\value{
A list containing the following elements:
  \item{eig}{a matrix containing all the eigenvalues, the percentage of variance and the cumulative percentage of variance}
  \item{var}{a list of matrices containing all the results for the categories (coordinates, square cosine, contributions, v.test)}
  \item{ind}{a list of matrices containing all the results for the products (coordinates, square cosine, contributions)}
  \item{group}{a list of matrices containing all the results for consumers (coordinates, square cosine, contributions)}
  \item{acm}{all the results of the MCA}
  \item{cooccur}{the reordered co-occurrence matrix among products}
  \item{reord}{the reordered matrix products*consumers}
  \item{cramer}{the Cramer's V matrix between all the consumers}
  \item{textual}{the results of the textual analysis for the products}
  \item{validity}{the elements of validity calculated for the first eigenvalue and the ellipses}
  \item{call}{a list with some statistics}
}

\references{
Cadoret, M., L\^e, S., Pag\`es, J. (2008) \emph{A novel Factorial Approach for analysing Sorting Task data}. 9th Sensometrics meeting. St Catharines, Canada\cr
Cadoret, M., L\^e, S., Pag\`es, J. (2009) \emph{A Factorial Approach for Sorting Task data (FAST)}. Food Quality and Preference.  20. pp. 410-417\cr
Cadoret, M., L\^e, S., Pag\`es, J. (2009) \emph{Missing values in categorization}. Applied Stochastic Models and Data Analysis (ASMDA). Vilnius, Lithuania\cr
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\examples{
\dontrun{
data(perfume)
## Example of FAST results
res.fast<-fast(perfume)
}
}
\keyword{multivariate}