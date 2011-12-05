\name{carto}
\alias{carto}

\title{Preference Mapping Techniques}

\description{
Performs preference mapping techniques based on multidimensional exploratory data analysis.
}

\usage{
carto(Mat, MatH, 
      level = 0, regmod = 1, coord = c(1, 2), asp = 1, 
      cex = 1.3, col = "steelblue4", font = 2, clabel = 0.8,
      label.j = FALSE, resolution = 200, nb.clusters = 0,
	    graph.tree=TRUE,graph.corr=TRUE,graph.carto=TRUE,
	    main=NULL,col.min=7.5,col.max=0)
}

\arguments{

\item{Mat}{a data frame corresponding to the axes of the map}
\item{MatH}{a data frame in which each row represent a product and each column represent the hedonic scores of a given 
consumer for the products}
\item{level}{the number of standard deviations used in the calculation of the preference response surface for all
the consumers}
\item{regmod}{the type of regression model used in the calculation of the preference response surface for all
the consumers. \code{regmod} = 1: quadratic model, 
\code{regmod} = 2: vector model, \code{regmod} = 3: circular model, \code{regmod} = 4: elliptical model}
\item{coord}{a vector of length 2, the rank of the axis used to display the results if \code{"manual"}
is not assigned to the \code{option} parameter}
\item{asp}{if 1 is assigned to that parameter, the graphic displays are output in an orthonormal coordinate system}
\item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
\item{col}{cf. function \code{\link{par}} in the \pkg{graphics} package}
\item{font}{cf. function \code{\link{par}} in the \pkg{graphics} package}
\item{clabel}{cf. the \pkg{ade4} package}
\item{label.j}{boolean, if T then the labels of the panelists who gave the hedonic scores are displayed}
\item{resolution}{resolution of the map}
\item{nb.clusters}{number of clusters to use (by default, 0 and the optimal numer of clusters is calculated}
\item{graph.tree}{boolean, if TRUE plots the tree in 2 dimensions}
\item{graph.corr}{boolean, if TRUE plots the variables factor map}
\item{graph.carto}{boolean, if TRUE plots the preference map}
\item{main}{an overall title for the plot}
\item{col.min}{define the color which match to the low levels of preference}
\item{col.max}{define the color which match to the high levels of preference}
}

\details{
The preference mapping methods are commonly used in the fields of market research 
and research and development to explore and understand the structure and tendencies of consumer preferences, 
to link consumer preference information to other data and to predict the behavior of consumers 
in terms of acceptance of a given product.\cr
This function refers to the method introduced by M. Danzart. A response surface is computed per consumer; then
according to certain threshold preference zones are delimited and finally superimposed.

}
\references{

Danzart M., Sieffermann J.M., Delarue J. (2004). New developments in preference mapping techniques: 
finding out a consumer optimal product, its sensory profile and the key sensory attributes. 
\emph{7th Sensometrics Conference, July 27-30, 2004, Davis, CA.}\cr

}

\author{
Francois Husson \email{husson@agrocampus-rennes.fr} \cr
Sebastien Le \email{Sebastien.Le@agrocampus-rennes.fr}
}

\seealso{\code{\link[FactoMineR]{MFA}}, \code{\link[FactoMineR]{GPA}}}

\examples{
## Example 1: carto for the sensory descriptors
data(cocktail)
res.pca <- PCA(senso.cocktail)
res.carto <- carto(res.pca$ind$coord[,1:2], hedo.cocktail)

## Example 2
\dontrun{
data(cocktail)
res.mfa <- MFA(cbind.data.frame(senso.cocktail,compo.cocktail),
    group=c(ncol(senso.cocktail),ncol(compo.cocktail)),
    name.group=c("senso","compo"))
res.carto <- carto(res.mfa$ind$coord[,1:2], hedo.cocktail)
}
}
\keyword{multivariate}
\keyword{models}
