\name{JAR}

\alias{JAR}

\title{JAR}

\description{Just About Right}

\usage{
JAR(x, col.p, col.j, col.pref, jarlevel="jar")
}

\arguments{
  \item{x}{data.frame}
  \item{col.p}{the position of the \emph{product} variable}
  \item{col.j}{the position of the \emph{panelist} variable}
  \item{col.pref}{the position of the \emph{preference} variable}
  \item{jarlevel}{a string corresponding to the jar level (the level must be the same for all the jar variables)}
}

\details{
Perform the penalty analysis. Two models are constructed.\cr
The one-dimensional model is constructed descriptor by descriptor. For descriptor_j the model is:\cr
Hedonic score = Descriptor_j_Not enough+ Descriptor_j_Too much \cr
The multi-dimensional model is constructed with all descriptors simultaneously:\cr
Hedonic score = Descriptor_1_Not enough+ Descriptor_1_Too much +...+ Descriptor_p_Not enough+ Descriptor_p_Too much+ Product + Judge\cr
}

\value{
Returns a list of 3 objects. \cr
The penalty1 object corresponds to the one-dimensional penalty results: a data-frame with the penalty coefficient in the first column, the standard deviation and the p-value for the test that the penalty is significantly different from 0.\cr
The penalty2 object corresponds to the mutli-dimensional penalty results: a data-frame with the penalty coefficient in the first column, the standard deviation and the p-value for the test that the penalty is significantly different from 0.
The Frequency object gives the percentage of times the non-jar categories are given for each product: a matrix with the non-jar categories in rows and the products in columns}

\seealso{ \code{\link{plot.JAR}}}

\examples{
\dontrun{
data(JAR)
res.jar <- JAR(JAR,col.p=13,col.j=1,col.pref=2)
plot(res.jar,name.prod="284", model=1)
 }
}

\author{Francois Husson}

\keyword{multivariate}
\keyword{dplot}
