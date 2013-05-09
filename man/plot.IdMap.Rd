\name{plot.IdMap}

\alias{plot.IdMap}

\title{Plot for the Ideal Mapping}

\description{
Option to replot by zooming on the ideal map.
}

\usage{
\method{plot}{IdMap}(x, xlim, ylim, levels.contour=NULL, color=FALSE, inverse=FALSE, \dots)
}

\arguments{
  \item{x}{An object of class \emph{IdMap}}
  \item{xlim}{The lower and upper limits of interest on the X-axis }
  \item{ylim}{The lower and upper limits of interest on the Y-axis}
  \item{levels.contour}{The levels (between 0 and 1) to consider for the colors on the surface plot.
  By default, they are set automatically based on the results }
  \item{color}{Boolean, define whether the map is in color or in black and white}
  \item{inverse}{Boolean, define whether the Black/While code should be inversed in the graph}
  \item{\dots}{further arguments passed to or from other methods}
}

\details{
Replot the Ideal Map by zooming on the area of interest
}

\value{
Return the IdMap zoomed on the area of interest.
}

\author{Worch Thierry (thierry@qistatistics.co.uk)}

\seealso{ \code{\link{IdMap}}}

\examples{
\dontrun{
data(perfume_ideal)
res <- IdMap(perfume_ideal, col.p=2, col.j=1, col.lik=ncol(perfume_ideal), id.recogn="id_")
plot.IdMap(res, xlim=c(-7,10), ylim=c(-5,7), levels.contour=NULL, color=TRUE)
plot.IdMap(res, xlim=c(-7,10), ylim=c(-5,7), levels.contour=NULL, color=FALSE, inverse=TRUE)
}
}

\keyword{dplot}
