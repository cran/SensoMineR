\name{senso.label}

\alias{senso.label}

\title{senso.label}

\description{
\code{scatter} is a generic function. It has methods for the classes \code{coa}, \code{dudi}, \code{fca}, \code{acm} and \code{pco}.\cr
The scale of the grid is situated on the right-top of the graph.\cr
The points are in the middle of the labels.\cr
This process plots the graphs of the multivariate analyses.\cr
The two axes have the same scale.\cr
}

\usage{
senso.label(dfxy, xax = 1, yax = 2, label = row.names(dfxy), clabel = 1, 
    pch = 20, cpoint = if (clabel == 0) 1 else 0, neig = NULL, 
    cneig = 2, xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE, 
    cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "", 
    csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL, 
    area = NULL, add.plot = FALSE) 
}
\arguments{
  \item{dfxy}{a data frame with two coordinates }
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{label}{a vector of strings of characters for the point labels}
  \item{clabel}{if not NULL, a character size for the labels, used with  \code{par("cex")*clabel}}
  \item{pch}{}
  \item{cpoint}{}
  \item{neig}{}
  \item{cneig}{}  
  \item{xlim}{}
  \item{ylim}{}
  \item{grid}{a logical value indicating whether a grid in the background of the plot should be drawn}
  \item{addaxes}{}
  \item{cgrid}{}
  \item{include.origin}{}
  \item{origin}{}
  \item{sub}{a string of characters to be inserted as legend}
  \item{csub}{a character size for the legend, used with \code{par("cex")*csub}}
  \item{possub}{a string of characters indicating the sub-title position ("topleft", "topright", "bottomleft", "bottomright")}
  \item{cgrid}{a character size, parameter used with par("cex")*\code{cgrid} to indicate the mesh of the grid}
  \item{pixmap}{}
  \item{contour}{}
  \item{area}{}
  \item{add.plot}{if TRUE uses the current graphics window}
}

\details{
The only difference between this function and the one programed by Daniel Chessel
is the frame displayed by the function \code{box}.
}

\author{Daniel Chessel \email{chessel@biomserv.univ-lyon1.fr}}

\keyword{internal}
