\name{plotpanelist}
\alias{plotpanelist}

\title{Plotpanelist}

\description{
Displays panelists' sensory profiles onto the products' space
}

\usage{
plotpanelist(mat, coord = c(1,2), name = FALSE, eig, cex = 1, color = NULL)
}

\arguments{
  \item{mat}{a data frame structured as the first element of the list resulting from the function construct.axes,
  i.e. the coordinates of the products with respect to the panel and to each panelists}
  \item{coord}{a length 2 vector specifying the components to plot}
  \item{name}{boolean, if T then the name of each panelist is displayed on the graph (by default, FALSE
  is assigned to that parameter)}
  \item{eig}{a matrix with the component of the factor analysis (in row) and the eigenvalues, the inertia and the cumulative inertia for each component. Typically, the \code{eig} output of the \code{construct.axes} function}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{color}{a vector with the colors used; by default there are 35 colors defined}
}

\value{
Returns a graph where each product is displayed with respect to a panel and to each panelist composing
the panel. Products described by the panel are displayed as square, they are displayed as circle when
they are described by each panelist.
}

\author{François Husson}

\examples{
data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
plotpanelist(axe$moyen, eig = signif(axe$eig,4))
}
\keyword{multivariate}
