\name{plot.WordCountAna}

\alias{plot.WordCountAna}

\title{Word-Count based methods Analysis (WordCountAna) graphs}

\description{
Draw the Word-Count based methods Analysis (WordCountAna) graphs
}

\usage{
\method{plot}{WordCountAna}(x,axes=c(1,2),choix="prod",lab=TRUE,color=NULL,pch=NULL,
    proba=0.05,xlim=NULL,ylim=NULL,cex=1,title=NULL,new.plot=TRUE,\dots)
}

\arguments{
  \item{x}{an object of class VerbAna}
  \item{axes}{a length 2 vector specifying the components to plot}
  \item{choix}{a string corresponding to the graph that you want to do ("prod" for the products graph, "panel" for the panellists graph, "dist" for the distinct-words graph, "cons" for the consensual words graph)}
  \item{lab}{boolean, if TRUE, the labels are drawn}
  \item{color}{the color to use to draw the graph}
  \item{pch}{either an integer specifying a symbol or a single character to be used as the default in plotting points. See points for possible values and their interpretation}
  \item{proba}{significance threshold considered to define consensual words (by default 0.05)}
  \item{xlim}{range for the plotted 'x' values, defaulting to the range of the finite values of 'x'}
  \item{ylim}{range for the plotted 'y' values, defaulting to the range of the finite values of 'y'}
  \item{cex}{numerical value giving the amount by which plotting text and symbols should be magnified relative to the default (by default 1)}
  \item{title}{string corresponding to the title of the graph you draw}
  \item{new.plot}{boolean, if TRUE, a new graphical device is created}
  \item{\dots}{further arguments passed to or from other methods, such as cex, cex.main, ...}
}

\value{
Returns one of the four following factor maps depending on the choice: products, panellists, distinct-words and consensual words.
}

\references{
Kostov, B., Becue-Bertaut, M., Husson, F., Pages, J., Cadoret, M., Torrens, J. and Urpi, P. (2012). A tool for detecting words with consensual meaning in verbalization tasks. 11th Sensometrics Conference, July 10-13, 2012, Rennes, France.
}

\author{Belchin Kostov \email{badriyan@clinic.ub.es}, Francois Husson \email{francois.husson@institut-agro.fr}, Monica Becue-Bertaut}

\seealso{ \code{\link{WordCountAna}}}

\examples{
data(perfume)
res<-WordCountAna(base=perfume,sep.word=";",graph=FALSE)
plot.WordCountAna(res,choix="prod")
plot.WordCountAna(res,choix="panel")
plot.WordCountAna(res,choix="dist")
plot.WordCountAna(res,choix="cons")
plot.WordCountAna(res,choix="cons",proba=0.1)
}

\keyword{multivariate}