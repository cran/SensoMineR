\name{search.desc.ktab}

\alias{search.desc.ktab}

\title{Search for discriminating descriptors}

\description{
This function is designed to select the significant descriptors in an object of class \code{ktab}.
}

\usage{
search.desc.ktab(ktab.donnee, level = 0.5)
}

\arguments{
  \item{ktab.donnee}{an object of class \code{ktab} (i.e. a list of data frames with the same row.names in common)}
  \item{level}{the threshold (P-value) below which variables are considered as discriminating for
  the analysis of variance model: \code{descriptor=Panelist+Product}}  
}

\details{
The input data set is an object of class \code{ktab} for which the first element of the list 
(i.e. the first data frame) must correspond to a 2
qualitative variables block, where the first variable refers to the \emph{panelist} variable
and the second to the \emph{product} variable.
}

\value{
Returns an object of class \code{ktab} with only discriminating variables.
}

\author{François Husson}

\examples{
data(chocolates)
#Construction of the ktab (J stands for Panelist and P stands for Product)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],
    chocolates[,4], chocolates[,5:ncol(chocolates)]), blocks = c(2,14), 
    tabnames = c("JP","Gr1"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee, level = 0.5)
}

\keyword{models}
