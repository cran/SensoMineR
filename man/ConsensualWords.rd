\name{ConsensualWords}

\alias{ConsensualWords}

\title{Consensual words for Sorting Task data}

\description{
This function is designed to point out the words that are used in a consensual way by consumers from a sorting task.
}
  
\usage{ConsensualWords(res.fast, nbtimes = 3, nbsimul = 500, proba = 0.05, 
    graph = TRUE, axes = c(1,2))}

\arguments{
  \item{res.fast}{an object of class fast}  
  \item{nbtimes}{minimum sample size for the word selection}
  \item{nbsimul}{the number of simulations used to compute Bootstrap}
  \item{proba}{the significance threshold considered to consider a word as consensual 
  (by default 0.05)}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{axes}{a length 2 vector specifying the components to plot}
}


\value{
A list containing the following elements:
  \item{Centroids}{coordinates of the consensual words on the dimensions of the fast result}
  \item{Within.inertia}{frequency of use of each word and within inertia associated}
  \item{Results.Bootstrap}{frequency of use of each word, within inertia associated and 
  p-value calculated according to the Bootstrap technique}
  \item{Consensual.words}{a list of significant consensual words sorted from the most 
  consensual to the less consensual}
}


\author{Francois Husson}

\examples{
\dontrun{
data(perfume)
## Example of FAST results
res.fast<-fast(perfume,sep.words=";")
res.consensual<-ConsensualWords(res.fast)
}
}
\keyword{multivariate}