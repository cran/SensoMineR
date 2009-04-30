\name{simulation}

\alias{simulation}

\title{Simulate virtual panels}

\description{Simulate virtual panels}

\usage{
simulation(axeAFM, nbchoix = NULL, nbgroup = 1, nbsimul = 500)
}

\arguments{
  \item{axeAFM}{matrix}
  \item{nbchoix}{Number of panelists in a virtual Panel}  
  \item{nbgroup}{Number of groups (by default 1 for a PCA)}  
  \item{nbsimul}{number of simulations}  
}

\details{
Called by the simulation procedure
}

\references{Husson, F., Le Dien, S. & Pagès, J.  (2005). Confidence ellipse for the sensory profiles obtained by Principal Components Analysis. Food Quality and Preference.  16 (3). pp. 245-250.
}

\author{François Husson}

\examples{
data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
simul <- simulation(axe)
plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))
}

\keyword{models}
\keyword{internal}
