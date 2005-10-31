\name{senso.corcircle}
\alias{senso.corcircle}

\title{Plot of the factorial maps of a correlation circle}

\description{
Performs the scatter diagram of a correlation circle. 
}
\usage{
senso.corcircle(dfxy, xax = 1, yax = 2, label = row.names(df),clabel = 1, grid = TRUE, sub = "", csub = 1, possub = "bottomleft",cgrid = 0, fullcircle = TRUE, box = FALSE, add.plot = FALSE,col="black",lty=1)
}
\arguments{
  \item{dfxy}{a data frame with two coordinates }
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{label}{a vector of strings of characters for the point labels}
  \item{clabel}{if not NULL, a character size for the labels, used with  \code{par("cex")*clabel}}
  \item{grid}{a logical value indicating whether a grid in the background of the plot should be drawn}
  \item{sub}{a string of characters to be inserted as legend}
  \item{csub}{a character size for the legend, used with \code{par("cex")*csub}}
  \item{possub}{a string of characters indicating the sub-title position ("topleft", "topright", "bottomleft", "bottomright")}
  \item{cgrid}{a character size, parameter used with par("cex")*\code{cgrid} to indicate the mesh of the grid}
  \item{fullcircle}{a logical value indicating whether the complete circle sould be drawn}
  \item{box}{a logical value indcating whether a box should be drawn}
  \item{add.plot}{if TRUE uses the current graphics window}
  \item{col}{color of the variables}
  \item{lty}{line type of the arrows}
}

\details{
The only difference between this function and the one programed by Daniel Chessel
is that \code{scatterutil.eti.circ} was replaced by \code{senso.scatterutil.eti.circ}
in order not to display frames around the names of the variables.
}

\author{Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}}

\keyword{multivariate}
\keyword{internal}
