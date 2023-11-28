\name{IdMapConsumer}

\alias{IdMapConsumer}

\title{Ideal Mapping by categories (IdMapConsumer)}

\description{
Create the ideal map and plot the ideal areas of the categories of qualitative variables. 
And perform 2 tests: a global test in order to highlight the significance of the difference between ideals of all the categories of the same variable; 
a pair comparison test to highlight the significance between 2 categories of the same variable.}

\usage{
IdMapConsumer(dataset.id, dataset.signa, col.p, col.j, col.lik, 
	num.col.var.signa, conf.level=0.95, id.recogn, nbchoix = NULL,
    nbsimul = 500, alpha = 0.05, coord = c(1, 2), precision = 0.1,
    levels.contour = NULL, color = FALSE, simusigni = 500) 
	}

\arguments{
  \item{dataset.id}{A matrix with at least two qualitative variables (\emph{consumer} and \emph{products})
  and a set of quantitative variables containing at least 2*A variables (for both \emph{perceived} and \emph{ideal} intensities)}
  \item{dataset.signa}{a data frame with n rows (individuals) and p columns (categorical variables)}
  \item{col.p}{The position of the \emph{product} variable in the dataframe \emph{dataset.id}}
  \item{col.j}{The position of the \emph{consumer} variable  in the dataframe \emph{dataset.id}}
  \item{col.lik}{The position of the \emph{liking} variable  in the dataframe \emph{dataset.id}}
  \item{id.recogn}{The sequence in the variable names which distinguish the ideal
  variables from the sensory variables. This sequence should be fixed and unique. \cr
  Each ideal variable should be preceeded by the corresponding perceived intensity variable.}
  \item{num.col.var.signa}{The position of the categorical variables in the dataframe \emph{dataset.signa} you want to plot the ideal area 
  of the different modalities/you want to know if the ideal product of the different modalities is significantly different}
  \item{conf.level}{Threshold used for the tests}
  \item{nbchoix}{The number of consumers forming a virtual panel, by default the number of panelists in the original panel}
  \item{nbsimul}{The number of simulations (corresponding to the number of virtual panels) used to compute the ellipses}
  \item{alpha}{The confidence level of the ellipses}
  \item{coord}{A length 2 vector specifying the components to plot}
  \item{precision}{The value defining the step when gridding the space}
  \item{levels.contour}{The levels (between 0 and 1) to consider for the colors on the surface plot.
  By default, they are set automatically based on the results}
  \item{color}{Boolean, define whether the map is in color or in black and white}
  \item{simusigni}{The number of simulations used to perform the global and the pair comparison test}
}

\details{
The IdMapConsumer, step by step:
Step 1: the classical IdMap is plotted with the method "ellipses"\cr
Step 2: for each modality of the categorical variable, the optimum of the ideal area is calculated with the method "density"\cr
Step 3: for each categorical variable given in \emph{num.col.var.signa}, simulations are performed giving the p-value 
for the global ant the pair comparison test. \cr
Step 4: if the global test is significant for a variable, the ideal areas of its modalities are plotted on the IdMap

This function needs the KernSmooth package.}

\value{
A list containing the following components:
  \item{PCA}{the results from the PCA used to create the sensory space}
  \item{idmap}{a list containing the results of the IdMap (\emph{data}), the weight for each consumer (\emph{j.weight})
  and the precision used.}
  \item{ideal}{a list containing the estimated profile of the ideal of reference (not available for the wIdMap) 
  as well as the percentage of consumers concerned}
  \item{coordobs}{The coordinates of all the ideals of all the categories on the sensory space}
  \item{test.global}{The results for the global test for each variables (observed inertia, critical inertia, P-value)}
  \item{test.paires}{The results for the pair comparison test  for each variables, between its ideal's categories(observed distance between two categories, critical distance, P-value)}

  The three last components are provided only if the user choose "color = FALSE", else no test and no ideal map with categories' ideal are performed.
 }

\references{Worch, T., Le, S., Punter, P., & Pages, J. (2012).
Construction of an Ideal Map (IdMap) based on the ideal profiles obtained directly from consumers.
\emph{Food Quality and Preference}, 26, 93-104.
}

\author{Melodie Sanchez, Sarah Sanchez, francois.husson@institut-agro.fr}

\seealso{ \code{\link{IdMap}}}

\examples{
\dontrun{
###Load the two datasets
data(cream_id)
data(cream_signa)

###Run the analysis and test the ideals of the variables from 1 to 12 
## for example with a confidence level of 90%, and 500 simulations.
res.idmap <- IdMapConsumer(cream_id, cream_signa, col.p=2, col.j=1, col.lik=29, 
   num.col.var.signa=c(1:12),conf.level=0.90,id.recogn="id_") 
}
}

\keyword{multivariate}