\name{hmfa}
\alias{hmfa}
\title{Hierarchical Multiple Factor Analysis}
\description{
Performs a hierarchical multiple factor analysis, using an object of class \code{list} of \code{data.frame}.
}
\usage{
hmfa(X, H, num = 6, coord = c(1,2))
}
\arguments{
  \item{X}{a \code{data.frame}}
  \item{H}{a list with one vector for each hierarchical level; in each vector the number of variables or the number of group
  constituting the group}
  \item{num}{number of graphs for the partial plots in one window}
  \item{coord}{a length 2 vector specifying the components to plot}
  }
\value{
  Returns a list including:
  \item{hmfa}{an object all the results of the hmfa}
  \item{coord.partial}{a list of lists of matrix}
}
\references{Le Dien, S. & Pagès, J. (2003) Hierarchical Multiple factor analysis: application to the comparison of sensory profiles, \emph{Food Quality and Preferences}, \bold{18 (6)}, 453-464.
}
\author{Sébastien Lê}

\details{

}

\examples{ 
data(sensopanels)
reshmfa<-hmfa(sensopanels, H = list(rep(14,7),c(4,3)))
}

\keyword{multivariate}
