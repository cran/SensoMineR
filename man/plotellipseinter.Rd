\name{plotellipseinter}
\alias{plotellipseinter}

\title{Plotellipseinter}
\description{Construit l'ensemble des ellipses de confiance}
\usage{
plotellipseinter(mat, alpha = 0.05, coord = c(1,2), nbgroup = 1,
    moy = TRUE, eig, cex = 1, color = NULL, title = NULL)
}

\arguments{
  \item{mat}{Matrix with all the points}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{alpha}{Confidence level of the ellipses}
  \item{nbgroup}{Confidence level of the ellipses}
  \item{moy}{Confidence level of the ellipses}
  \item{eig}{a matrix with the component of the factor analysis (in row) and the eigenvalues, the inertia and the cumulative inertia for each component. Typically, the \code{eig} output of the \code{construct.axes} function}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
}

\author{François Husson}

\keyword{internal}
