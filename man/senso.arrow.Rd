\name{senso.arrow}

\alias{senso.arrow}

\title{senso.arrow}

\description{
\code{scatter} is a generic function. It has methods for the classes \code{coa}, \code{dudi}, \code{fca}, \code{acm} and \code{pco}.\cr
The scale of the grid is situated on the right-top of the graph.\cr
The points are in the middle of the labels.\cr
This process plots the graphs of the multivariate analyses.\cr
The two axes have the same scale.\cr
}

\usage{
senso.arrow(dfxy, xax = 1, yax = 2, label =row.names(dfxy), 
    clabel = 1, pch = 20, cpoint = 0, edge = TRUE, 
    origin = c(0, 0), xlim = NULL, ylim = NULL, grid = TRUE, 
    addaxes = TRUE, cgrid = 1, sub = "", csub = 1.25, 
    possub = "bottomleft", pixmap = NULL, contour = NULL,
    area = NULL, add.plot = FALSE,col="black",lty=1)
}
\arguments{

  \item{dfxy}{a data frame with two coordinates }
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{label}{a vector of strings of characters for the point labels}
  \item{clabel}{if not NULL, a character size for the labels, used with  \code{par("cex")*clabel}}
  \item{pch}{if \code{cpoint > 0}, an integer specifying the symbol or the single character to be used in plotting points}
  \item{cpoint}{a character size for plotting the points, used with par("cex")*cpoint. If zero, no points are drawn}
  \item{edge}{a logical value indicating whether the arrows should be plotted}
  \item{origin}{the fixed point in the graph space, by default c(0,0) the origin of axes. The arrows begin at cent}  
  \item{xlim}{the ranges to be encompassed by the x-axis, if NULL they are computed}
  \item{ylim}{the ranges to be encompassed by the y-axis, if NULL they are computed}
  \item{grid}{a logical value indicating whether a grid in the background of the plot should be drawn}
  \item{addaxes}{a logical value indicating whether the axes should be plotted}
  \item{cgrid}{a character size, parameter used with par("cex")*cgrid, to indicate the mesh of the grid}
  \item{sub}{a string of characters to be inserted as legend}
  \item{csub}{a character size for the legend, used with \code{par("cex")*csub}}
  \item{possub}{a string of characters indicating the sub-title position ("topleft", "topright", "bottomleft", "bottomright")}
  \item{pixmap}{an object 'pixmap' displayed in the map background}
  \item{contour}{a data frame with 4 columns to plot the contour of the map : each row gives a segment (x1,y1,x2,y2)}
  \item{area}{a data frame of class 'area' to plot a set of surface units in contour}
  \item{add.plot}{if TRUE uses the current graphics window}
  \item{col}{color of the arrows}
  \item{lty}{line type of the arrows}
}

\author{Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}}

\keyword{internal}
