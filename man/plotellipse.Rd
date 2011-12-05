\name{plotellipse}
\alias{plotellipse}

\title{Plot confidence ellipses}

\description{Plot confidence ellipses.}

\usage{
plotellipse(mat, alpha = 0.05, coord = c(1,2), eig, cex = 1, 
    color = NULL, title = NULL)
}

\arguments{
  \item{mat}{Matrix with all the points}
  \item{alpha}{the confidence level of the ellipses}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{eig}{a matrix with the component of the factor analysis (in row) and the eigenvalues, the inertia and the cumulative inertia for each component. Typically, the \code{eig} output of the \code{construct.axes} function}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
}
\author{François Husson}

\examples{
\dontrun{
data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
simul <- simulation(axe)
plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))
#######################################
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, group = c(6,8), 
    name.group = c("A-F","T-S"),scale.unit = TRUE)
simul <- simulation(axe, nbgroup = (ncol(axe$partiel)-2)/(ncol(axe$moyen)-2))
plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))
}}
\keyword{dplot}
\keyword{internal}
