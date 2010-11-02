\name{plot.fast}

\alias{plot.fast}

\title{Make Factorial Approach for Sorting Task data (FAST) graphs}

\description{
Plot the graphs for Factorial Approach for Sorting Task data (FAST).
}

\usage{
plot.fast(x,choix="ind", axes = c(1, 2), xlim = NULL, ylim = NULL, invisible = NULL,
    col.ind = "blue", col.var = "red", col.quali.sup = "darkred", 
    col.ind.sup ="darkblue", col.quanti.sup = "black",label = "all",
    cex = 1,lab.grpe = TRUE, title = NULL, habillage = "none", palette = NULL,
    new.plot = TRUE, \dots)
    }

\arguments{
  \item{x}{an object of class fast}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{choix}{the graph to plot ("ind" for the products, "var" for the vcategories, "group" for the consumers)}
  \item{xlim}{range for the plotted 'x' values, defaulting to the range of the finite values of 'x'}
  \item{ylim}{range for the plotted 'y' values, defaulting to the range of the finite values of 'y'}
  \item{habillage}{give no color for the individuals ("none"), or color the products among a consumer (give the number of the consumer)}
  \item{col.ind}{a color for the products}
  \item{col.var}{a color for the categories}
  \item{col.quali.sup}{a color for the supplementary categories}
  \item{col.ind.sup}{a color for the supplementary individuals}
  \item{col.quanti.sup}{a color for the quantitative supplementary variables}
  \item{label}{a list of character for the elements which are labelled (by default, all the elements are labelled ("ind", "var"))}
  \item{invisible}{string indicating if some points should not be drawn ("ind" or "var")}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{lab.grpe}{boolean, if TRUE, the consumers are labelled}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{palette}{the color palette used to draw the points. By default colors are chosen. If you want to define the colors : palette=palette(c("black","red","blue")); or you can use: palette=palette(rainbow(30)), or in black and white for example: palette=palette(gray(seq(0,.9,len=25)))}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
Returns the products factor map, the categories factor map and the consumers factor map.
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fast}}}

\examples{
\dontrun{
data(perfume)
res.fast <- fast(perfume,graph=FALSE)
plot.fast(res.fast,choix="ind",invisible="var",habillage=5)
plot.fast(res.fast,choix="group")
}
}

\keyword{dplot}
