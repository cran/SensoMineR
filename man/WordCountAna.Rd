\name{WordCountAna}

\alias{WordCountAna}

\title{Word-Count based methods Analysis  (WordCountAna)}

\description{
Sensory methods as labelled sorting task, check-all-that-apply (CATA), ultra-flash profiling (UFP) and open-ended questions can be used to collect free-text descriptions of products through word-count based methods. A data frame with rows-products and columns-panellists is considered for the analysis. WordCountAna performs a multiple factor analysis for contingency tables keeping all the information in the comparison of the products. The identification of the consensual words which have the same meaning for most of the panellists eases the interpretation of the word-count based methods and solves the problems arising from the large diversity of vocabulary as the different meanings possibly associated to a same word. A test, based on resampling techniques, allows for assessing the significance of the consensus.
}

\usage{
WordCountAna (base, sep.word = NULL, ncp = Inf, nb.panel = 3, nb.simul = 500, 
    proba = 0.05, graph = TRUE, axes = c(1,2))
}

\arguments{
  \item{base}{a data frame with n rows (products) and p columns (panellists). Each cell corresponds to a free-text description used to describe a product by a panellist}
  \item{sep.word}{a string with all the characters which correspond to separator of words (by default, NULL and is considered equal to "; (),?./:'!$=+;{}<>[]@-")}
  \item{ncp}{number of dimensions kept in the results and to compute the within-inertia}
  \item{nb.panel}{minimum number of panellists who used the same word in order to define consensual words (by default 3)}
  \item{nb.simul}{number of bootstrap simulations (by default 500)}
  \item{proba}{significance threshold considered to define consensual words (by default 0.05)}
  \item{graph}{boolean, if TRUE a graph is displayed}
  \item{axes}{a length 2 vector specifying the components to plot}
}

\value{
  \item{mfact}{a list of matrices containing all the results for multiple factor analysis for contingency tables}
  \item{dist.words}{a matrix containing the results for distinct words (number of times that used and number of panellists that pronounced)}
  \item{centroids}{a matrix containing the coordinates of the centroids of distinct-words}
  \item{cons}{a matrix containing the results of bootstrap resampling for distinct-words pronounced by at least "nb.panel" panellists (number of times that used, number of panellists that pronounced and the significance of the consensus)}
  \item{cons.words}{a vector of consensual words assessed by bootstrap resampling}

Returns the products factor map, panellists factor map, distinct-words factor map and consensual words factor map.
}

\references{
Becue-Bertaut, M. and Pages, J. (2008). Multiple factor analysis and clustering of a mixture of quantitative, categorical and frequency data. \emph{Computational Statistice and Data Analysis}, 52, 3255-3268.
Kostov, B., Becue-Bertaut, M. and Husson, F. (2012). Multiple Factor Analysis for Contingency Tables in FactoMineR Package. \emph{The R journal}
Kostov, B., Becue-Bertaut, M., Husson, F., Pages, J., Cadoret, M., Torrens, J. and Urpi, P. (2012). A tool for detecting words with consensual meaning in verbalization tasks. 11th Sensometrics Conference, July 10-13, 2012, Rennes, France.
}

\author{Belchin Kostov \email{badriyan@clinic.ub.es}, Francois Husson \email{husson@agrocampus-ouest.fr}, Monica Becue-Bertaut}

\seealso{ \code{\link{textual}}, \code{\link{MFA}}, \code{\link{plot.WordCountAna}}}

\examples{
data(perfume)
res<-WordCountAna(base=perfume,sep.word=";")
}

\keyword{multivariate}