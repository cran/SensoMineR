\name{graphinter}

\alias{graphinter}

\title{Graphical display of the interaction between two qualitative variables}

\description{
This function is designed to display the interaction between two qualitative variables, in most cases
the \emph{product} and the \emph{session} variables.
}

\usage{
graphinter(donnee, col.p, col.j, firstvar, lastvar=ncol(donnee), 
    numr = 2,numc = 2)
}

\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}  
  \item{col.p}{the position of one categorical variables of interest (the \emph{product} variable)}
  \item{col.j}{the position of one categorical variables of interest (the \emph{session} variable)}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee})}
  \item{numr}{the number of graphs per row (by default 2)}
  \item{numc}{the number of graphs per column (by default 2)}
}


\details{
The data set must be balanced (or not unbalanced too much).
}

\value{
If the variables of interest are the \emph{product} and the \emph{session} variables, 
a list containing the following components:
  \item{prod}{a data frame of dimension (\emph{p,q}), the means over the panelists and the sessions for the
  \emph{p} products and the \emph{q} sensory descriptors}
  \item{seance}{as many matrices of dimension (\emph{p,q}) as there are sessions, 
  the means over the panelists for the \emph{p} products, the \emph{q} sensory descriptors and for each session}
The graphical display of the interaction for each sensory descriptor.
}

\references{P. Lea, T. Naes, M. Rodbotten. \emph{Analysis of variance for sensory data}. \cr
H. Sahai, M. I. Ageel. \emph{The analysis of variance}.
}

\author{
François Husson, Sébastien Lê}

\seealso{ \code{\link{aov}}}

\examples{
\dontrun{
data(chocolates)
graphinter(sensochoc, col.p = 4, col.j = 2, firstvar = 5, lastvar = 12,
    numr = 1, numc = 1)
}
}
\keyword{models}
