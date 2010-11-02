\name{cpa}
\alias{cpa}

\title{Consumers' Preferences Analysis}

\description{
Performs preference mapping techniques based on multidimensional exploratory data analysis. This methodology is
oriented towards consumers' preferences; here consumers are pictured according only to their preferences. 
In this manner, the distance between two consumers is very natural and easy to interpret, and a clustering of 
the consumers is also very easy to obtain. 
}

\usage{
cpa(senso, hedo, coord=c(1,2), center = TRUE, scale = TRUE, 
    nb.clusters = 0, scale.unit = FALSE, 
    col = terrain.colors(45)[1:41])
}

\arguments{
\item{senso}{a data frame of dimension (\emph{p,k}), where \emph{p} is the number of products and \emph{k} the
number of sensory descriptors}
\item{hedo}{a data frame of dimension (\emph{p,j}), where \emph{p} is the number of products and \emph{j} the
number of consumers or panelists}
\item{coord}{a length 2 vector specifying the components to plot}
\item{center}{boolean, if TRUE then data are mean centered}
\item{scale}{boolean, if TRUE then data are scaled to unit variance}
\item{nb.clusters}{number of clusters to use (by default, 0 and the optimal numer of clusters is calculated}
\item{scale.unit}{boolean, if TRUE then PCA is made on scaled data}
\item{col}{color palette}
}

\details{
This methodology is oriented towards consumers' preferences; here, consumers are pictured according only to their preferences. 
In this manner, the distance between two consumers is very natural and easy to interpret, 
and a clustering of the consumers is also very easy to obtain using a classic hierarchical 
clustering procedure performed on Euclidian distances with the Ward's minimum variance criterion. 
The originality of the representation is that the characteristics of the products 
are also superimposed to the former picture.
}

\value{

Return the following results:
\item{clusters}{the cluster number allocated to each consumer}
\item{result}{the coordinates of the panelists, of the clusters, of the archetypes}
\item{prod.clusters}{a list with as many elements as there are clusters; each element of the list gathers
the specific products for its corresponding cluster}
\item{desc.clusters}{the correlation coefficients between the average hedonic scores per cluster and the sensory descriptors}
A dendogram which highlight the clustering, a correlation circle that displays the hedonic scores, a graph of the consumers such as two consumers are all 
the more close that they do like the same products, as many graphs as there are variables: for a given variable, 
each consumer is colored according to the coefficient of correlation based on his hedonic scores and the variable.
}

\references{
S. Lê, F. Husson, J. Pagès (2005). Another look at sensory data: how to "have your salmon and eat it, too!". 
\emph{6th Pangborn sensory science symposium, August 7-11, 2005, Harrogate, UK.}
}

\author{
François Husson \email{François.Husson@agrocampus-rennes.fr} \cr
Sébastien Lê \email{Sebastien.Le@agrocampus-rennes.fr}
}

\examples{
\dontrun{
data(cocktail)
res.cpa = cpa(cbind(compo.cocktail, senso.cocktail), hedo.cocktail)
## If you prefer a graph in black and white and with 3 clusters
res.cpa = cpa(cbind(compo.cocktail, senso.cocktail), hedo.cocktail, 
    name.panelist = TRUE, col = gray((50:1)/50), nb.clusters = 3)
}
}
\keyword{multivariate}
