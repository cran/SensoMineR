\name{IdMap}

\alias{IdMap}

\title{Ideal Mapping (IdMap)}

\description{
Create the ideal map, a map based on the ideal profiles provided by the consumers.
}

\usage{
IdMap(dataset, col.p, col.j, col.lik, id.recogn, nbchoix=NULL, nbsimul=500,
    alpha=0.05, coord=c(1,2), precision=0.1, levels.contour=NULL, color=FALSE, 
	cons.eq=FALSE)
}

\arguments{
  \item{dataset}{A matrix with at least two qualitative variables (\emph{consumer} and \emph{products})
  and a set of quantitative variables containing at least 2*A variables (for both \emph{perceived} and \emph{ideal} intensities)}
  \item{col.p}{The position of the \emph{product} variable }
  \item{col.j}{The position of the \emph{consumer} variable }
  \item{col.lik}{The position of the \emph{liking} variable }
  \item{id.recogn}{The sequence in the variable names which distinguish the ideal
  variables from the sensory variables. This sequence should be fixed and unique. \cr
  Each ideal variable should be preceeded by the corresponding perceived intensity variable. }
  \item{nbchoix}{The number of consumers forming a virtual panel, by default the number of panelists in the original panel}
  \item{nbsimul}{The number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{alpha}{The confidence level of the ellipses}
  \item{coord}{A length 2 vector specifying the components to plot}
  \item{precision}{The value defining the step when gridding the space}
  \item{levels.contour}{The levels (between 0 and 1) to consider for the colors on the surface plot.
  By default, they are set automatically based on the results}
  \item{color}{Boolean, define whether the map is in color or in black and white}
  \item{cons.eq}{Boolean, define whether the IdMap (by default) or the wIdMap is performed}
}

\details{
  The IdMap, step by step:\cr
  Step 1: the sensory and ideal variables are separated into two tables.\cr
  Step 2: the product space is created by PCA on the averaged sensory table (averaged by product).\cr
  Step 3: the averaged ideal product of each consumer is projected as supplementary entities in this space.\cr
  Step 4: confidence ellipses are created around each individual averaged ideal product using truncated total bootstrap.\cr
  Step 5: for each consumer, the space is grid and the position where the ideal area
  is defined is marked: individual surfaces of response are created.\cr
  Step 6: (optional) the ellipses can be balanced by applying individual weight
  (all the ellipses have a weigth of 1, however the size of the ellipse). wIdMap is then performed.\cr
  Step 7: all the individual surface plots are added together and a surface plot is created.\cr
}

\value{
A list containing the following components:
  \item{PCA}{the results from tthe PCA used to create the sensory space}
  \item{idmap}{a list containing the results of the IdMap (\emph{data}), the weight for each consumer (\emph{j.weight})
  and the precision used.}
  \item{ideal}{a list caining the estimtated profile of the ideal of reference (not available for the wIdMap) 
  as well as the percentage of consumers concerned}
}

\references{Worch, T., Lê, S., Punter, P., & Pagès, J. (2012).
Construction of an Ideal Map (IdMap) based on the ideal profiles obtained directly from consumers.
\emph{Food Quality and Preference}, 26, 93-104.
}

\author{Worch Thierry (thierry@qistatistics.co.uk)}

\seealso{ \code{\link{plot.IdMap}}, \code{\link{carto}}, \code{\link{boot}}}

\examples{
\dontrun{
data(perfume_ideal)
#! For the IdMap
res.IdMap <- IdMap(perfume_ideal, col.p=2, col.j=1, col.lik=ncol(perfume_ideal), 
     id.recogn="id_", cons.eq=FALSE)
plot.IdMap(res.IdMap, x=c(-7,10), y=c(-5,7), levels.contour=NULL, color=TRUE)

#! For the wIdMap
res.wIdMap <- IdMap(perfume_ideal, col.p=2, col.j=1, col.lik=ncol(perfume_ideal), 
     id.recogn="id_", cons.eq=TRUE)
}
}

\keyword{multivariate}