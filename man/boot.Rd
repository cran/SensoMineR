\name{boot}

\alias{boot}

\title{Simulate virtual panels for several functions}

\description{Simulate virtual panels for the sorting task, the napping, the sorting napping, 
the free choice profiling, the hierarchical sorting task}

\usage{
boot(X, method = "sorting", axes = 1:2, scale = TRUE, ncp = NULL, group = NULL,
    nbsim = 200,level.conf = 0.95,nbchoix = NULL,color = NULL,cex = 0.8, 
    title = NULL, new.plot = TRUE)
}

\arguments{
  \item{X}{data.frame}
  \item{method}{String with the method to use. The argument can be "sorting" (the default for sorting task data), 
"napping" (for napping data), "sortnapping" (for sorted napping), 
"freechoice" (for free choice profiling), "hsort" (for hierarchical sorting task data).}  
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{scale}{boolean, used when method="freechoice"; if TRUE, the variables are scaled}
  \item{ncp}{number of components used to procrustes the virtual subspaces on the true subspace; by default NULL and the number of components is estimated}
  \item{group}{a list indicating the number of variables in each group; used when method="freechoice" or method="hsort"}
  \item{nbsim}{the number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{level.conf}{confidence level used to construct the ellipses. By default, 0.95}
  \item{nbchoix}{the number of panelists forming a virtual panel, by default the number of panelists in the original panel}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
}

\details{
Calculate virtual panels by bootstrap of the panelists. For each virtual panel, calculate the mean 
configuration and procrustes this configuration on the true configuration obtained from the true panel.
}

\value{
Returns a list with estim.ncp which corresponds to the output of the estim_ncp function (function which estimates the number of components) and the simul object which can be used with the plotellipse function.
}

\author{Marine Cadoret and François Husson}

\examples{
######## Napping example
data(napping)
res <- boot(napping.don,method="napping")

######## Sorting task example
data(perfume)
res <- boot(perfume,method="sorting")

######## Sorted task napping example
data(smoothies)
res <- boot(smoothies,method="sortnapping")

######## Hierarchical sorting task example
data(cards)
group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
res <- boot(cards,method="hsort", group=group.cards)

######## Free choice profiling example
data(perfume_fcp)
res <- boot(perfume_fcp, method="freechoice", group = c(12,7,7,7,6,8))
}

\keyword{multivariate}
\keyword{dplot}
