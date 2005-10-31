\name{mfasenso}

\alias{mfasenso}

\title{Make MFA with illustrative individuals}

\description{

}

\usage{
mfasenso(ktab,ktab.illu=NULL,scale.unit=TRUE,nbcoord=2,poids=NULL)
}

\arguments{
  \item{ktab}{ktableau des individus actifs}
  \item{ktab.illu}{ktableau des individus illustratifs}  
  \item{scale.unit}{Booléen qui prend la valeur vraie si les variables sont normées}  
  \item{nbcoord}{number of coordinates}  
  \item{poids}{pour mettre des poids aux variables (intérêt limité très certainement)}  
}

\details{
Programme d'analyse factorielle: 2 ktableaux le 1er actif, le 2ème illustratif
ressort les coord des points moyens et partiels, des individus actifs et illustratifs
Si dans les k-tableaux il n'y a qu'un groupe, fait aussi des ACP.
}

\value{

}

\references{
Escofier B., Pagès J. (1988, 1990, 1993, 1998) Analyses factorielles simples et multiples. Objectifs méthodes et interprétation. Dunod. Paris.
Pagès, J. & Husson, F. (2005) Multiple Factor Analysis with confidence ellipses: a useful methodology to study the relations between sensory and instrumental data. Journal of chemometrics.
}

\author{François Husson}

\examples{
data(chocolates)
ktab.donnee <- ktab.data.frame(cbind.data.frame(chocolates[,1],chocolates[,4],chocolates[,-(1:4)]),blocks=c(2,14),tabnames=c("JP","Gr1"))
}

\keyword{dplot}
\keyword{internal}
