\name{plotellipse}
\alias{plotellipse}

\title{Plot confidence ellipses}

\description{Plot confidence ellipses.}

\usage{
plotellipse(mat,alpha=0.05,coord=c(1,2))
}

\arguments{
  \item{mat}{Matrix with all the points}
  \item{alpha}{the confidence level of the ellipses}
  \item{coord}{a length 2 vector specifying the components to plot}
}


\author{François Husson}


\examples{
\dontrun{
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,14),tabnames=c("JP","Gr1"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=0.5)
axe <- construct.axes(ktab.interesting.desc,scale.unit=TRUE)
simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbsimul=500)
plotellipse (simul,alpha=0.05)
#######################################
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,6,4,4),tabnames=c("JP","A","F","T"))
ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=0.5)
axe <- construct.axes(ktab.interesting.desc,scale.unit=TRUE)
simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbsimul=500)
plotellipse (simul,alpha=0.05)
}
}
\keyword{dplot}
\keyword{internal}
