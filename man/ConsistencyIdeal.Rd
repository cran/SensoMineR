\name{ConsistencyIdeal}

\alias{ConsistencyIdeal}

\title{Sensory and Hedonic consistency of the ideal data}

\description{
Evaluate the sensory and hedonic consistency of the ideal data, both at the consumer and panel level.
}

\usage{

ConsistencyIdeal(dataset, col.p, col.j, col.lik, id.recogn, type="both", scale.unit=TRUE,
                  ncp=NULL, axes=c(1,2), nbsim=0, replace.na=FALSE, graph=TRUE)
}

\arguments{
  \item{dataset}{A matrix with at least two qualitative variables (\emph{consumer} and \emph{products})
  and a set of quantitative variables containing at least 2*A variables (for both \emph{perceived} and \emph{ideal} intensities)}
  \item{col.p}{The position of the \emph{product} variable }
  \item{col.j}{The position of the \emph{consumer} variable }
  \item{col.lik}{The position of the \emph{liking} variable }
  \item{id.recogn}{The sequence in the variable names which distinguish the ideal
  variables from the sensory variables. This sequence should be fixed and unique. \cr
  Each ideal variable should be preceeded by the corresponding perceived intensity variable. }
  \item{type}{Define whether you want the sensory consistency only ("sensory"), the hedonic consistency only ("hedonic"), or both ("both")}
  \item{scale.unit}{Boolean, if TRUE the descriptors are scaled to unit variance}
  \item{ncp}{Number of dimensions kept in the results}
  \item{axes}{A length 2 vector specifying the components to plot}
  \item{nbsim}{The number of simulations performed. By default (=0), no simulations are performed and only the results for the real data are given}
  \item{replace.na}{Boolean, define whether the missing values (in the correlation matrix calculated for the consistency at the consumer level) should be ignored or replaced by 0}
  \item{graph}{Boolean, define whether the distribution of the correlation coefficient should be plot}
}

\details{
  SENSORY CONSISTENCY\cr
  A the panel level:\cr
  A PCA is performed on the table crossing the J consumers in rows and the A ideal variables
  in columns (the averaged or corrected averaged is then considered).\cr
  On this space, the sensory description of the P products (P rows) on the A attributes
  is projected as supplementary entities while the hedonic table crossing
  the J consumers (in rows) and the P products (in columns) is projected as supplementary variables.\cr
  The sensory consistency is measured by the correspondence between the same products seen
  through the sensory and through the hedonic descriptions.\cr
  At the consumer level:\cr
  For each consumer, the correlation between the (corrected) ideal ratings and
  the correlation between the hedonic scores and the perceived intensity of each attribute is calculated.\cr
  A test on this correlation coefficient is performed for each consumer.\cr
  The distribution of these correlations coefficients are also given graphically.
  \cr
  \cr
  HEDONIC CONSISTENCY \cr
  For each consumer, a PCR-model expressing the liking scores in function of the perceived intensity is created.\cr
  Once the model is created, the model is applied to the ideal ratings provided by the consumer
  considered and the hedonic score of the ideal product is estimated.\cr
  This hedonic score is then compared to the hedonic scores provided to the products tested.\cr
  If simulations are asked, the same procedure is estimated after re-sampling the vector of hedonic scores.\cr
  In that case, the distribution of the estimated ideal hedonic score can be estimated
  under H0 and the significance of the estimated ideal hedonic score can be done.\cr
}

\value{
A list containing the results for the sensory and hedonic consistency:\cr
  \item{Senso}{contains the results of the sensory consistency}
  \item{Senso$panel}{results for the consistency at the panel level including:}
  \item{Senso$panel$dataset}{the datasets used for the different PCA}
  \item{Senso$panel$PCA.ideal}{the results of the PCA for the creation of the ideal space}
  \item{Senso$panel$PCA.ideal_hedo}{the results of the PCA with projection of the hedonic scores}
  \item{Senso$panel$PCA.ideal_senso}{the results of the PCA with the proojection of the sensory descriptions}
  \item{Senso$panel$correlation}{the correlation between the \emph{product} projected from the \emph{sensory} and \emph{hedonic} points of view}
  \item{Senso$conso}{results of the consistency at the consumer level including: }
  \item{Senso$conso$driver.lik}{the linear drivers of liking (correlation between perceived intensity and liking score for each attribute)}
  \item{Senso$conso$correlations}{the correlations between drivers of liking and the difference (ideal-perceived) intensity}
\cr
  \item{Hedo$R2}{the R2 coefficients of the indivvidual models}
  \item{Hedo$hedo}{a list containing the hedonic scores for the \emph{product}, \emph{ideal products}, \emph{average ideal product}
  and the \emph{standardized ideal product} for each consumer}
  \item{Hedo$simulation}{a list including the estimated hedonic score for each simulation, the p-value and the matrix of simulations used}
}

\references{
{Worch, T., Le, S., Punter, P., & Pages, J. (2012).
Assessment of the consistency of ideal profiles according to non-ideal data for IPM.
\emph{Food Quality and Preference}, 24, 99-110.},
{Worch, T., Le, S., Punter, P., & Pages, J. (2012).
Extension of the consistency of the data obtained with the Ideal Profile Method:
Would the ideal products be more liked than the tested products?
\emph{Food Quality and Preference}, 26, 74-80.}
}

\author{Thierry Worch (thierry@qistatistics.co.uk)}

\seealso{ \code{\link{panelperf}}, \code{\link{paneliperf}}}

\examples{
\dontrun{
data(perfume_ideal)
res <- ConsistencyIdeal(perfume_ideal, col.p=2, col.j=1, 
   col.lik=ncol(perfume_ideal), id.recogn="id_", 
   type="both", nbsim=100)
}
}
\keyword{models}