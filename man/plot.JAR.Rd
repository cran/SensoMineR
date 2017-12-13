\name{plot.JAR}

\alias{plot.JAR}

\title{Make penalty graph}

\description{
Plot the graphs for the penalty analysis.
}

\usage{
\method{plot}{JAR}(x, name.prod, model=1, confidence=TRUE, level=0.05, \dots)
    }

\arguments{
  \item{x}{an object of class fast}
  \item{name.prod}{a string corresonding to the name of the product for which penalty is drawn}
  \item{model}{an integer: 1 means that the one-dimensional penalties are drawn, 2 means that the multidimensional penalties are drawn}
  \item{confidence}{A boolean: if true, the confidence interval is drawn for each penalty}
  \item{level}{threshold for the selection of the categories that are plotted}
  \item{\dots}{further arguments passed to or from other methods, such as cex, cex.main, ...}
}

\value{
Draw a plot of the penalty against the ferquency of the categories.
}

\author{Francois Husson}

\seealso{ \code{\link{JAR}}}

\examples{
\dontrun{
data(JAR)
res.jar <- JAR(JAR,col.p=13,col.j=1,col.pref=2)
plot.JAR(res.jar,name.prod="284", model=1)
 }
}

\keyword{dplot}
