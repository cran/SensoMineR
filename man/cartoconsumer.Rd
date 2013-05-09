\name{cartoconsumer}
\alias{cartoconsumer}

\title{Preference Mapping Techniques and segmentation of consumers}

\description{
Performs preference mapping techniques based on multidimensional exploratory data analysis and segmentation of consumers.
}

\usage{
cartoconsumer(res, data.pref, nb.clust=0, seuil=0.8, consol=TRUE, ncp=5,
		 scale.conso=TRUE,graph.carto=TRUE,graph.hcpc=FALSE, graph.group=FALSE,
		col.min=7.5, col.max=0, contrast=0.2, level=0, asp=0,lwd=2)
				}

\arguments{

\item{res}{the result of a factor analysis}
\item{data.pref}{a data frame in which each row represent a product and each column represent the hedonic scores of a given 
consumer for the products}
\item{nb.clust}{	an integer. If 0, the tree is cut at the level the user clicks on. If -1, the tree is automatically cut at the suggested level (see details). If a (positive) integer, the tree is cut with nb.cluters clusters}
\item{seuil}{the size of the area kept for each group of consumers}
\item{consol}{a boolean. If TRUE, a k-means consolidation is performed}
\item{ncp}{number of dimensions kept in the results (by default 5)}
\item{scale.conso}{scale data by consumer}
\item{graph.carto}{if TRUE, the preference map is displayed. If FALSE, no graph is displayed}
\item{graph.hcpc}{if TRUE, graphics of segmentation (trees and indivuals map) are displayed. If FALSE, no graph are displayed}
\item{graph.group}{if TRUE, preference maps for each group are displayed. If FALSE, no map are displayed}
\item{col.min}{define the color which match to the low levels of preference}
\item{col.max}{define the color which match to the high levels of preference}
\item{contrast}{define the color contrast between groups' areas and the rest of the map}
\item{level}{the number of standard deviations used in the calculation of the preference response surface for all
the consumers}
\item{asp}{if 1 is assigned to that parameter, the graphic displays are output in an orthonormal coordinate system}
\item{lwd}{control the line width for the outlines defining groups' areas}
}

\details{
The preference mapping methods are commonly used in the fields of market research 
and research and development to explore and understand the structure and tendencies of consumer preferences, 
to link consumer preference information to other data and to predict the behavior of consumers 
in terms of acceptance of a given product.\cr
This function refers to the method introduced by M. Danzart. A segmentation of consumers is performed, and a 
preference map is displayed for each group of consumers. The original preference map is built, the areas of 
each group are underlined thanks to a contrast, and the number of consumers is shown.

}
\references{
Danzart M., Sieffermann J.M., Delarue J. (2004). New developments in preference mapping techniques: 
finding out a consumer optimal product, its sensory profile and the key sensory attributes. 
\emph{7th Sensometrics Conference, July 27-30, 2004, Davis, CA.}\cr

}

\author{
Francois Husson \email{husson@agrocampus-rennes.fr} \cr
Sophie Birot and Celia Pontet
}

\seealso{\code{\link[FactoMineR]{MFA}}, \code{\link[FactoMineR]{GPA}}, \code{\link[SensoMineR]{carto}}}

\examples{
\dontrun{
## Example 1: carto on the sensory descriptors
data(cocktail)
res.pca <- PCA(senso.cocktail)
results1 <- cartoconsumer(res.pca, hedo.cocktail)
results2 <- cartoconsumer(res.pca, hedo.cocktail,
      graph.hcpc=TRUE,graph.group=TRUE)
}

## Example 2
\dontrun{
data(cocktail)
res.mfa <- MFA(cbind.data.frame(senso.cocktail,compo.cocktail),
    group=c(ncol(senso.cocktail),ncol(compo.cocktail)),
    name.group=c("senso","compo"))
results3 <- cartoconsumer(res.mfa, hedo.cocktail)
}
}
\keyword{multivariate}
\keyword{models}
\keyword{segmentation}