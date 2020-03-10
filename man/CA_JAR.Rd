\name{CA_JAR}

\alias{CA_JAR}

\title{Make a correspondence analysis on the JAR data}

\description{
Plot the CA graph.
}

\usage{
CA_JAR(x, col.p, col.j, col.pref, jarlevel="jar")
}

\arguments{
  \item{x}{data.frame}
  \item{col.p}{the position of the \emph{product} variable}
  \item{col.j}{the position of the \emph{panelist} variable}
  \item{col.pref}{the position of the \emph{preference} variable}
  \item{jarlevel}{a string corresponding to the jar level (the level must be the same for all the jar variables)}
}

\value{
Draw a CA graph with the preference data as supplementary qualitative variables, the products as rows, and the categories of the jar variables as columns
}

\author{Francois Husson}

\examples{
\dontrun{
data(JAR)
res <- CA_JAR(x = JAR, col.p = 13, col.j = 1, col.pref = 2, jarlevel = "jar")
plot.CA(res$res.CA, invisible = "row", cex = 0.8)
 }
}

\keyword{dplot}
