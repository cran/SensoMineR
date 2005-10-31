\name{gpa}

\alias{gpa}

\title{Generalised Procrustes Analysis}

\description{
Performs a Generalised Procrustes Analysis (GPA) that takes into account missing values.
}
\usage{
gpa(df,tolerance=10^-10,nbiteration=200,scal=TRUE,coord=c(1,2))
}

\arguments{
  \item{df}{an object of class \code{ktab} (i.e. a list of data frames with the same row.names in common)}
  \item{tolerance}{a threshold with respect to which the algorithm stops, i.e. when the difference between 
   the GPA loss function at step \emph{n} and \emph{n+1} is less than \code{tolerance}}
  \item{nbiteration}{the maximum number of iterations until the algorithm stops}
  \item{scal}{a boolean, if TRUE (which is the default value) scaling is required}
  \item{coord}{a length 2 vector specifying the components to plot}
}

\details{
Performs a Generalised Procrustes Analysis (GPA) that takes into account missing values:
some data frames of \code{df} may have non described or non evaluated rows, i.e. rows with missing 
values only.\cr
The algorithm used here is the one developed by Commandeur.
}


\value{
A list containing the following components:
\item{RV}{a matrix of RV coefficients between partial configurations}
\item{RVs}{a matrix of standardized RV coefficients between partial configurations}
\item{simi}{a matrix of Procrustes similarity indexes between partial configurations}
\item{scaling}{a vector of isotropic scaling factors}
\item{consensus}{a matrix of consensus configuration}
\item{Xfin}{an array of partial configurations}
\item{PANOVA}{a list of "Procrustes Analysis of Variance" tables, per assesor (config), per product(objet), per dimension (dimension)}
}

\references{
Commandeur, J.J.F (1991) \emph{Matching configurations}.DSWO press, Leiden University.\cr
Dijksterhuis, G. & Punter, P. (1990) Interpreting generalized procrustes analysis "Analysis of Variance " tables, 
\emph{Food Quality and Preference}, \bold{2}, 255--265\cr
Gower, J.C (1975) Generalized Procrustes analysis, \emph{Psychometrika}, \bold{40}, 33--50 \cr
Kazi-Aoual, F., Hitier, S., Sabatier, R., Lebreton, J.-D., (1995) Refined approximations to permutations tests 
for multivariate inference. Computational Statistics and Data Analysis, \bold{20}, 643--656 \cr
Qannari, E.M., MacFie, H.J.H, Courcoux, P. (1999) Performance indices
and isotropic scaling factors in sensory profiling, \emph{Food Quality and Preference}, \bold{10}, 17--21
}

\author{Elisabeth Morand}

\examples{
data(sensopanels)
ktabsenso <- ktab.data.frame(sensopanels, blocks = rep(14,7))
gpa(ktabsenso)
}
  
\keyword{multivariate}
