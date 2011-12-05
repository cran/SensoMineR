\name{plot.fahst}

\alias{plot.fahst}

\title{Make Factorial Approach for Hierarchical Sorting Task data (FAHST) graphs}

\description{
Plot the graphs for Factorial Approach for Hierarchical Sorting Task data (FAHST).
}

\usage{
plot.fahst(x,choix="ind", axes = c(1, 2), xlim = NULL, ylim = NULL, invisible = NULL,
    col.ind = "blue", col.var = "red", lab.ind=TRUE,lab.var=TRUE, cex = 1,
    lab.lev=TRUE,lab.grpe = TRUE, title = NULL, habillage = "none", habillage.lev = "none",
    traj = FALSE, palette = NULL, new.plot = TRUE, \dots)
    }

\arguments{
  \item{x}{an object of class fahst}
  \item{choix}{the graph to plot ("ind" for the products and the categories, "group" for the consumers and "level" for the levels)}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{xlim}{range for the plotted 'x' values, defaulting to the range of the finite values of 'x'}
  \item{ylim}{range for the plotted 'y' values, defaulting to the range of the finite values of 'y'}
  \item{invisible}{string indicating if some points should not be drawn ("ind" or "var")}
  \item{col.ind}{a color for the products}
  \item{col.var}{a color for the categories}
  \item{lab.ind}{boolean, if TRUE, the products are labelled}
  \item{lab.var}{boolean, if TRUE, the categories associated with categorization are labelled}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{lab.lev}{boolean, if TRUE, the levels are labelled}
  \item{lab.grpe}{boolean, if TRUE, the consumers are labelled}
  \item{title}{string corresponding to the title of the graph you draw (by default NULL and a title is chosen)}
  \item{habillage}{give no color for the individuals ("none"), or color the products according to one of the levels of a consumer (give the number of the colomn corresponding to the level)}
  \item{habillage.lev}{give no color for the levels ("none"), color the levels according to consumer ("subject") or color the levels according to the number of the level ("level")}  
  \item{traj}{boolean, if TRUE, trajectories are drawn between levels of the same consumer}
  \item{palette}{the color palette used to draw the points. By default colors are chosen. If you want to define the colors : palette=palette(c("black","red","blue")); or you can use: palette=palette(rainbow(30)), or in black and white for example: palette=palette(gray(seq(0,.9,len=25)))}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{\dots}{further arguments passed to or from other methods}
}

\value{
Returns the products factor map, the categories factor map, the levels factor map and the consumers factor map.
}

\author{Marine Cadoret, S\'ebastien L\^e \email{sebastien.le@agrocampus-ouest.fr}}

\seealso{ \code{\link{fahst}}}

\examples{
\dontrun{
data(cards)
group<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
res.fahst <- fahst(cards,group,graph=FALSE)
plot.fahst(res.fahst,choix="ind",invisible="var",habillage=2,title="Cards colored according to level 2 of subject 1")
plot.fahst(res.fahst,choix="level",traj=TRUE)
}
}
\keyword{dplot}
