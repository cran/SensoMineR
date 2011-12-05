\name{fcp}

\alias{fcp}

\title{Free choice profiling}

\description{Free choice profiling with confidence ellipses}

\usage{
fcp(X,group,scale=TRUE, ncp = NULL, axes=c(1,2), name.group = NULL, level.conf = 0.95, nbsim=500, nbchoix=NULL, cex=1, color=NULL, title=NULL, new.plot=TRUE, graph=c("ind","var","ellipse"))
}

\arguments{
  \item{X}{data.frame}
  \item{group}{a list indicating the number of variables in each group; used when method="freechoice" or method="hsort"}
  \item{scale}{boolean, used when method="freechoice"; if TRUE, the variables are scaled}
  \item{ncp}{number of components used to procrustes the virtual subspaces on the true subspace; NULL by default and the number of components is estimated}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{name.group}{the names of each group of variables}
  \item{level.conf}{confidence level used to construct the ellipses. By default, 0.95}
  \item{nbsim}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{nbchoix}{the number of panelists forming a virtual panel, by default the number of panelists in the original panel}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{graph}{list with the graphs to draw; "ind", "var" and "ellipse" by default}
}

\details{
Perform MFA on the data frame and calculate confidence ellipses around the products
}

\value{
Returns a list with the result of the MFA and the bootstraped results that can be used with the plotellipse function.
}

\examples{
\dontrun{
data(perfume_fcp)
res <- fcp(perfume_fcp, group = c(12,7,7,7,6,8))
}
}

\author{François Husson}

\keyword{multivariate}
\keyword{dplot}
