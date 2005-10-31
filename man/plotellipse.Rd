\name{plotellipse}
\alias{plotellipse}

\title{Plot confidence ellipses}

\description{Plot confidence ellipses.}

\usage{
plotellipse(mat, alpha = 0.05, coord = c(1,2), eig, cex = 1, color = NULL)
}

\arguments{
  \item{mat}{Matrix with all the points}
  \item{alpha}{the confidence level of the ellipses}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{eig}{a matrix with the component of the factor analysis (in row) and the eigenvalues, the inertia and the cumulative inertia for each component. Typically, the \code{eig} output of the \code{construct.axes} function}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
}
\author{François Husson}

\examples{
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1], 
    chocolates[,4], chocolates[,-(1:4)]), blocks = c(2,14), 
    tabnames = c("JP","Gr1"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee, level = 0.5)
axe <- construct.axes(ktab.interesting.desc, scale.unit = TRUE)
simul <- simulation(axe, nbbloc = length(ktab.interesting.desc$blo)-1)
plotellipse (simul, alpha = 0.05, eig = axe$eig)
#######################################
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],
    chocolates[,4],chocolates[,-(1:4)]), blocks = c(2,6,4,4), 
    tabnames = c("JP","A","F","T"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee, level = 0.5)
axe <- construct.axes(ktab.interesting.desc, scale.unit = TRUE)
simul <- simulation(axe, nbbloc = length(ktab.interesting.desc$blo)-1)
plotellipse (simul, alpha = 0.05, eig = axe$eig)
}
\keyword{dplot}
\keyword{internal}
