\name{coltable}

\alias{coltable}

\title{Color the cells of a data frame according to 2 threshold levels}

\description{

Return a colored display of a data frame according to 2 threshold levels.

}

\usage{
coltable(matrice,nbrow=nrow(matrice),nbcol=ncol(matrice),level.lower=0.05,col.lower="mistyrose",level.upper=1e10,col.upper="lightblue",cex=0,nbdec=4,main.title=NULL)
}

\arguments{
  \item{matrice}{a data frame (or a matrix) with only quantitative variables}
  \item{nbrow}{the number of rows to be displayed (by default, \code{nrow(matrice)})}
  \item{nbcol}{the number of columns to be displayed (by default, \code{ncol(matrice)})}
  \item{level.lower}{the threshold below which cells are colored in \code{col.lower}}
  \item{col.lower}{the color used for \code{level.lower}}
  \item{level.upper}{the threshold above which cells are colored in \code{col.upper}}
  \item{col.upper}{the color used for \code{level.upper}}
  \item{cex}{cf. function \code{\link{par}} in the \pkg{graphics} package}
  \item{nbdec}{the number of decimal places displayed}
  \item{main.title}{title of the graph(s)}
}

\details{
This function is very useful especially when there are a lot of values to check.
}

\author{
François Husson, Sébastien Lê
}

\examples{
## Example 1
data(chocolates)
resdecat<-decat(chocolates,formul="~Product+Panelist",firstvar=5,graph=FALSE)
coltable(magicsort(resdecat$tabT),level.lower=-1.96,level.upper=1.96,main.title="Products' description")

## Example 2
data(chocolates)
resaverage<-averagetable(chocolates,formul="~Product+Panelist",firstvar=5)
coltable(magicsort(resaverage),level.upper=6,level.lower=4,main.title="Average by chocolate")

## Example 3
data(chocolates)
resperf<-paneliperf(chocolates,formul="~Product+Panelist+Product:Panelist",formul.j="~Product",col.j=1,firstvar=5,lastvar=12,synthesis=FALSE,graph=FALSE)
resperfprob<-magicsort(resperf$prob.ind,method="median")
coltable(resperfprob,level.lower=0.05,level.upper=1,main.title="P-value of the F-test (by panelist)")

resperfr2<-magicsort(resperf$r2.ind,method="median",ascending=FALSE)
coltable(resperfr2,level.lower=0.00,level.upper=0.85,main.title="Adjusted R-square (by panelist)")
}

\keyword{color}
