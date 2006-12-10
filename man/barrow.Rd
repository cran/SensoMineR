\name{barrow}

\alias{barrow}

\title{Barplot per row with respect to a set of quantitative variables}

\description{
Returns as many barplots as there are rows in a matrix. The barplots are automatically generated 
for all the quantitative variables.
}

\usage{
barrow(donnee, numr = 2, numc = 2, numchar = 8, color = "lightblue",
      title = NULL)
}

\arguments{
  \item{donnee}{a data frame of dimension (\emph{p,q}), where p is the number of products
  and \emph{q} is the number of sensory descriptors for instance}
  \item{numr}{the number of barplots to be displayed per row (by default 2)}
  \item{numc}{the number of barplots to be displayed per column (by default 2)}
  \item{numchar}{the number of character used to write the boxplot labels (by default 8)}
  \item{color}{the color of the barplots (by default "lightblue")}
  \item{title}{the title used in the graphs}
}

\details{
Missing values are ignored when forming barplots. 
}

\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.\cr
Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey, P. A. (1983) \emph{Graphical Methods for Data Analysis}. Wadsworth & Brooks/Cole. 
}

\author{ 
Sébastien Lê \email{Sebastien.Le@agrocampus-rennes.fr}
}

\seealso{ 
\code{\link{plot}}
}

\examples{
data(chocolates)
resdecat<-decat(sensochoc, formul = "~Product+Panelist", firstvar = 5, 
    graph = FALSE)
\dontrun{
barrow(resdecat$tabT)
barrow(resdecat$coeff, color = "orange")
}
}
\keyword{univar}
