\name{paneliperf}

\alias{paneliperf}

\title{Panelists' performance according to their capabilities to dicriminate between products}

\description{
Computes automatically P-values, Vtests, residuals, r-square for 
each category of a given qualitative variable (e.g. the \emph{panelist} variable);\cr
Computes he agreement between each panelist and the panel results;\cr
Gives the panel results (optional).
}

\usage{
paneliperf(donnee, formul, formul.j = "~Product", col.j, firstvar,
      lastvar = ncol(donnee), synthesis = FALSE, random = TRUE, 
      graph = FALSE)
}

\arguments{
  \item{donnee}{a data frame made up of at least two qualitative variables 
  (\emph{product}, \emph{panelist}) and a set of quantitative variables (sensory descriptors)}
  \item{formul}{the aov model used for the panel}
  \item{formul.j}{the aov model used for each panelist (no \emph{panelist} effect allowed)}
  \item{col.j}{the position of the \emph{panelist} variable}
  \item{firstvar}{the position of the first endogenous variable}
  \item{lastvar}{the position of the last endogenous variable (by default the last column of \code{donnee}}
  \item{synthesis}{boolean, the possibility to have the anova results for the panel model}
  \item{random}{boolean, the status of the Panelist variable in the anova model for the panel}
  \item{graph}{boolean, draws the PCA and MFA graphs}
}

\details{
  The \code{formul} parameter must be filled in by an analysis of variance model and must 
  begin with the categorical variable of interest (e.g. the product effect)
  followed by the different other factors of interest (and their combinations).
  E.g.:\code{formul = "~Product+Session"}.
}

\value{
A list containing the following components:
  \item{prob.ind}{a matrix which rows are the panelist, which columns are the endogenous variables (in most cases the sensory descriptors) and
  which entries are the P-values associated to the AOV model}
  \item{vtest.ind}{a matrix which rows are the panelist, which columns are the endogenous variables (in most cases the sensory descriptors) and
  which entries are the Vtests associated to the AOV model}
  \item{res.ind}{a matrix which rows are the panelist, which columns are the endogenous variables (in most cases the sensory descriptors) and
  which entries are the residuals associated to the AOV model}
    \item{r2.ind}{a matrix which rows are the panelist, which columns are the endogenous variables (in most cases the sensory descriptors) and
  which entries are the R-square associated to the AOV model}
    \item{signif.ind}{a vector with the number of significant descriptors per panelist}
    \item{agree.ind}{a matrix with as many rows as there are panelists and as many columns as there are descriptors and 
    the entries of this matrix are the correlation coefficients between the product coefficients for the panel 
   and for the panelists}
\cr
  \item{complete}{a matrix with the v-test corresponding to the p.value (see \code{p.values} below), the median of the agreement (see \code{agree} upper), the standard deviation of the panel anova model (see \code{res} below)}
  \item{p.value}{a matrix of dimension (\emph{k,m}) of \emph{P-values} associated with the F-test for the 
  \emph{k} descriptors and the \emph{m} factors and their combinations considered in the analysis
  of variance model of interest}
  \item{variability}{a matrix of dimension (\emph{k,m}) where the entries correspond to the percentages of variability
  due to the effects introduced in the analysis of variance model of interest}
  \item{res}{a vector of dimension \emph{k} of residual terms for the analysis
  of variance model of interest}
  \item{r2}{a vector of dimension \emph{k} of r-squared for the analysis of variance model of interest}

The usual graphs when MFA is performed on the data.frame
resulting from vtest.ind and agree.ind. \cr
The PCA graphs for the complete output.\cr
}

\references{P. Lea, T. Naes, M. Rodbotten. \emph{Analysis of variance for sensory data}.
H. Sahai, M. I. Ageel. \emph{The analysis of variance}.
}

\author{François Husson, Sébastien Lê}

\seealso{ \code{\link{panelperf}}, \code{\link{aov}}}

\examples{
\dontrun{
data(chocolates)
res<-paneliperf(sensochoc, formul = "~Product+Panelist+Session+
  Product:Panelist+Product:Session+Panelist:Session",
  formul.j = "~Product", col.j = 1, firstvar = 5, synthesis = TRUE)
resprob<-magicsort(res$prob.ind, method = "median")
coltable(resprob, level.lower = 0.05, level.upper = 1,
    main.title = "P-value of the F-test (by panelist)")
hist(resprob,main="Histogram of the P-values",xlab="P-values")

resr2<-magicsort(res$r2.ind, method = "median", ascending = FALSE)
coltable(resr2, level.lower = 0.00, level.upper = 0.85,
    main.title = "Adjusted R-square (by panelist)")

resagree<-magicsort(res$agree, sort.mat = res$r2.ind, method = "median")
coltable(resagree, level.lower = 0.00, level.upper = 0.85,
    main.title = "Agreement between panelists")
hist(resagree,main="Histogram of the agreement between panelist and panel",
    xlab="Correlation coefficient between the product effect for 
    panelist and panel")

coltable(magicsort(res$p.value, sort.mat = res$p.value[,1], bycol = FALSE,
    method = "median"),
    main.title = "Panel performance (sorted by product P-value)")
}
}

\keyword{models}
