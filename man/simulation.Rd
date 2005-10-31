\name{simulation}

\alias{simulation}

\title{Simulate virtual panels}

\description{

}

\usage{
simulation(axeAFM,nbchoix=NULL,nbbloc,nbsimul=500)
}

\arguments{
  \item{axeAFM}{matrix}
  \item{nbchoix}{Number of panelist in a virtual Panel}  
  \item{nbbloc}{Number of bloc for the ktable}  
  \item{nbsimul}{number of simulations}  
}

\details{
Called by the simulation procedure
}

\value{

}

\references{Husson, F., Le Dien, S. & Pagès, J.  (2005). Confidence ellipse for the sensory profiles obtained by Principal Components Analysis. Food Quality and Preference.  16 (3). pp. 245-250.
}

\author{François Husson}

\examples{
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,14),tabnames=c("JP","Gr1"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=0.5)
axe <- construct.axes(ktab.interesting.desc,scale.unit=TRUE)
simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbsimul=500)
}

\keyword{models}
\keyword{internal}
