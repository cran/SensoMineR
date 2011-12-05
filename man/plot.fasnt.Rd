\name{plot.fasnt}

\alias{plot.fasnt}

\title{Make Factorial Approach for Sorting Napping Task data (FASNT) graphs}

\description{
Plot the graphs for Factorial Approach for Sorting Napping Task data (FASNT).
}

\usage{
plot.fasnt(x,choix="ind", axes = c(1, 2), xlim = NULL, ylim = NULL, invisible = NULL,
    col.ind = "blue", col.var = "red", lab.ind=TRUE,lab.var=TRUE, lab.coord=TRUE, lab.partial=TRUE,
    cex = 1,lab.grpe = TRUE, title = NULL, habillage = "none", palette = NULL,
    new.plot = TRUE, \dots)
    }




\arguments{
  \item{x}{an object of class fast}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{choix}{the graph to plot ("ind" for the products and the categories, "group" for the consumers and "partial" for the partial products)}
  \item{xlim}{range for the plotted 'x' values, defaulting to the range of the finite values of 'x'}
  \item{ylim}{range for the plotted 'y' values, defaulting to the range of the finite values of 'y'}
  \item{habillage}{give no color for the individuals ("none"), or color the products among a consumer (give the number of the consumer)}
  \item{col.ind}{a color for the products}
  \item{col.var}{a color for the categories}
  \item{lab.ind}{boolean, if TRUE, the products are labelled}
  \item{lab.var}{boolean, if TRUE, the categories associated with categorization are labelled}
  \item{lab.coord}{boolean, if TRUE, the napping variables are labelled}
  \item{lab.partial}{boolean, if TRUE, the partial products are labelled}
  \item{invisible}{string indicating if some points should not be drawn ("ind" or "var")}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{lab.grpe}{boolean, if TRUE, the consumers are labelled}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{palette}{the color palette used to draw the points. By default colors are chosen. If you want to define the colors : palette=palette(c("black","red","blue")); or you can use: palette=palette(rainbow(30)), or in black and white for example: palette=palette(gray(seq(0,.9,len=25)))}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
Returns the products factor map, the categories factor map, the coordinates factor map and the consumers factor map.
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fasnt}}}

\examples{
data(smoothies)
res.fasnt <- fasnt(smoothies, first="nappe",graph=FALSE)
plot.fasnt(res.fasnt,choix="ind",invisible="var",habillage=15,
  title="Objects colored according to the groups provided by consumer 5")
plot.fasnt(res.fasnt,choix="partial",lab.partial=FALSE)
}

\keyword{dplot}
