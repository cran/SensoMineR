pkgname <- "SensoMineR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('SensoMineR')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ardi")
### * ardi

flush(stderr()); flush(stdout())

### Name: ardi
### Title: Automatic Research of DIvergences between scores
### Aliases: ardi
### Keywords: univar

### ** Examples

## Not run: 
##D data(chocolates)
##D ardi(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
## End(Not run)
  


cleanEx()
nameEx("averagetable")
### * averagetable

flush(stderr()); flush(stdout())

### Name: averagetable
### Title: Computes a (products,descriptors) matrix
### Aliases: averagetable
### Keywords: models

### ** Examples

data(chocolates)
resaverage<-averagetable(sensochoc, formul = "~Product+Panelist",
    firstvar = 5)
coltable(magicsort(resaverage), level.upper = 6,level.lower = 4,
    main.title = "Average by chocolate")

res.pca = PCA(resaverage, scale.unit = TRUE)



cleanEx()
nameEx("barrow")
### * barrow

flush(stderr()); flush(stdout())

### Name: barrow
### Title: Barplot per row with respect to a set of quantitative variables
### Aliases: barrow
### Keywords: univar

### ** Examples

data(chocolates)
resdecat<-decat(sensochoc, formul = "~Product+Panelist", firstvar = 5, 
    graph = FALSE)
## Not run: 
##D barrow(resdecat$tabT)
##D barrow(resdecat$coeff, color = "orange")
## End(Not run)



cleanEx()
nameEx("boot")
### * boot

flush(stderr()); flush(stdout())

### Name: boot
### Title: Simulate virtual panels for several functions
### Aliases: boot
### Keywords: multivariate dplot

### ** Examples

## Not run: 
##D ######## Napping example
##D data(napping)
##D res <- boot(napping.don,method="napping")
##D 
##D ######## Sorting task example
##D data(perfume)
##D res <- boot(perfume,method="sorting")
##D 
##D ######## Sorted task napping example
##D data(smoothies)
##D res <- boot(smoothies,method="sortnapping")
##D 
##D ######## Hierarchical sorting task example
##D data(cards)
##D group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D res <- boot(cards,method="hsort", group=group.cards)
##D 
##D ######## Free choice profiling example
##D data(perfume_fcp)
##D res <- boot(perfume_fcp, method="freechoice", group = c(12,7,7,7,6,8))
## End(Not run)



cleanEx()
nameEx("boxprod")
### * boxprod

flush(stderr()); flush(stdout())

### Name: boxprod
### Title: Boxplot per category with respect to a categorical variable and
###   a set of quantitative variables
### Aliases: boxprod
### Keywords: univar

### ** Examples

data(chocolates)
boxprod(sensochoc, col.p = 4, firstvar = 5, numr = 2, numc = 2)



cleanEx()
nameEx("cards")
### * cards

flush(stderr()); flush(stdout())

### Name: cards
### Title: Cards
### Aliases: cards
### Keywords: datasets

### ** Examples

## Not run: 
##D data(cards)
##D ## Example of FAHST
##D group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D res.fahst<-fahst(cards,group=group.cards)
## End(Not run)



cleanEx()
nameEx("carto")
### * carto

flush(stderr()); flush(stdout())

### Name: carto
### Title: Preference Mapping Techniques
### Aliases: carto
### Keywords: multivariate models

### ** Examples

## Example 1: carto for the sensory descriptors
data(cocktail)
res.pca <- PCA(senso.cocktail)
res.carto <- carto(res.pca$ind$coord[,1:2], hedo.cocktail)

## Example 2
## Not run: 
##D data(cocktail)
##D res.mfa <- MFA(cbind.data.frame(senso.cocktail,compo.cocktail),
##D     group=c(ncol(senso.cocktail),ncol(compo.cocktail)),
##D     name.group=c("senso","compo"))
##D res.carto <- carto(res.mfa$ind$coord[,1:2], hedo.cocktail)
## End(Not run)



cleanEx()
nameEx("cartoconsumer")
### * cartoconsumer

flush(stderr()); flush(stdout())

### Name: cartoconsumer
### Title: Preference Mapping Techniques and segmentation of consumers
### Aliases: cartoconsumer
### Keywords: multivariate models segmentation

### ** Examples

## Not run: 
##D ## Example 1: carto on the sensory descriptors
##D data(cocktail)
##D res.pca <- PCA(senso.cocktail)
##D results1 <- cartoconsumer(res.pca, hedo.cocktail)
##D results2 <- cartoconsumer(res.pca, hedo.cocktail,
##D       graph.hcpc=TRUE,graph.group=TRUE)
## End(Not run)

## Example 2
## Not run: 
##D data(cocktail)
##D res.mfa <- MFA(cbind.data.frame(senso.cocktail,compo.cocktail),
##D     group=c(ncol(senso.cocktail),ncol(compo.cocktail)),
##D     name.group=c("senso","compo"))
##D results3 <- cartoconsumer(res.mfa, hedo.cocktail)
## End(Not run)



cleanEx()
nameEx("chocolates")
### * chocolates

flush(stderr()); flush(stdout())

### Name: chocolates
### Title: Chocolates data
### Aliases: chocolates
### Keywords: datasets

### ** Examples

data(chocolates)
decat(sensochoc, formul = "~Product+Panelist", firstvar = 5, graph = FALSE)



cleanEx()
nameEx("cocktail")
### * cocktail

flush(stderr()); flush(stdout())

### Name: cocktail
### Title: Cocktail data
### Aliases: cocktail
### Keywords: datasets

### ** Examples

data(cocktail)



cleanEx()
nameEx("coltable")
### * coltable

flush(stderr()); flush(stdout())

### Name: coltable
### Title: Color the cells of a data frame according to 4 threshold levels
### Aliases: coltable
### Keywords: color

### ** Examples

## Example 1
data(chocolates)
resdecat<-decat(sensochoc, formul = "~Product+Panelist", firstvar = 5,
    graph = FALSE)
resaverage<-averagetable(sensochoc, formul = "~Product+Panelist", 
    firstvar = 5)
resaverage.sort = resaverage[rownames(magicsort(resdecat$tabT)),
    colnames(magicsort(resdecat$tabT))]
coltable(resaverage.sort, magicsort(resdecat$tabT), 
    level.lower = -1.96, level.upper = 1.96,
    main.title = "Average by chocolate")

## Example 3
## Not run: 
##D data(chocolates)
##D resperf<-paneliperf(sensochoc, 
##D     formul = "~Product+Panelist+Product:Panelist", 
##D     formul.j = "~Product", col.j = 1, firstvar = 5, lastvar = 12, 
##D     synthesis = FALSE, graph = FALSE)
##D resperfprob<-magicsort(resperf$prob.ind, method = "median")
##D coltable(resperfprob, level.lower = 0.05, level.upper = 1, 
##D     main.title = "P-value of the F-test (by panelist)")
##D 
##D resperfr2<-magicsort(resperf$r2.ind, method = "median", 
##D     ascending = FALSE)
##D coltable(resperfr2, level.lower = 0.00, level.upper = 0.85, 
##D     main.title = "Adjusted R-square (by panelist)")
## End(Not run)



cleanEx()
nameEx("compo.cocktail")
### * compo.cocktail

flush(stderr()); flush(stdout())

### Name: compo.cocktail
### Title: Composition of the cocktails data
### Aliases: compo.cocktail
### Keywords: datasets

### ** Examples

data(cocktail)



cleanEx()
nameEx("construct.axes")
### * construct.axes

flush(stderr()); flush(stdout())

### Name: construct.axes
### Title: Coordinates of individuals and illustrative individuals for PCA
###   or MFA
### Aliases: construct.axes
### Keywords: multivariate

### ** Examples

## Example1: PCA 
data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
 
## Example2: MFA (two groups of variables)
data(chocolates) 
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, group = c(6,8), 
    name.group = c("A-F","T-S"),scale.unit = TRUE)



cleanEx()
nameEx("cpa")
### * cpa

flush(stderr()); flush(stdout())

### Name: cpa
### Title: Consumers' Preferences Analysis
### Aliases: cpa
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(cocktail)
##D res.cpa = cpa(cbind(compo.cocktail, senso.cocktail), hedo.cocktail)
##D ## If you prefer a graph in black and white and with 3 clusters
##D res.cpa = cpa(cbind(compo.cocktail, senso.cocktail), hedo.cocktail, 
##D     name.panelist = TRUE, col = gray((50:1)/50), nb.clusters = 3)
## End(Not run)



cleanEx()
nameEx("decat")
### * decat

flush(stderr()); flush(stdout())

### Name: decat
### Title: DEscription of CATegories
### Aliases: decat
### Keywords: models

### ** Examples

### Example 1
data(chocolates)
## model (AOV): " descriptor = product + panelist "
resdecat<-decat(sensochoc, formul="~Product+Panelist", firstvar = 5)
barrow(resdecat$tabT)
barrow(t(resdecat$tabT), numr = 3, numc = 3)
barrow(resdecat$coeff, color = "orange") 

### Example 2
data(chocolates)
## model (AOV): " descriptor = product + panelist "
res2 <-decat(sensochoc, formul="~Product+Panelist", firstvar = 5,
    proba=1, graph = FALSE)



cleanEx()
nameEx("fahst")
### * fahst

flush(stderr()); flush(stdout())

### Name: fahst
### Title: Factorial Approach for Hierarchical Sorting Task data
### Aliases: fahst
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(cards)
##D ## Example of FAHST results
##D group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D res.fahst<-fahst(cards,group=group.cards)
## End(Not run)



cleanEx()
nameEx("fasnt")
### * fasnt

flush(stderr()); flush(stdout())

### Name: fasnt
### Title: Factorial Approach for Sorting Napping Task data
### Aliases: fasnt
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(smoothies)
##D ## Example of FASNT results
##D res.fasnt<-fasnt(smoothies,first="nappe",sep.word=";")
## End(Not run)



cleanEx()
nameEx("fast")
### * fast

flush(stderr()); flush(stdout())

### Name: fast
### Title: Factorial Approach for Sorting Task data
### Aliases: fast
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(perfume)
##D ## Example of FAST results
##D res.fast<-fast(perfume,sep.words=";")
## End(Not run)



cleanEx()
nameEx("fcp")
### * fcp

flush(stderr()); flush(stdout())

### Name: fcp
### Title: Free choice profiling
### Aliases: fcp
### Keywords: multivariate dplot

### ** Examples

## Not run: 
##D data(perfume_fcp)
##D res <- fcp(perfume_fcp, group = c(12,7,7,7,6,8))
## End(Not run)



cleanEx()
nameEx("graphinter")
### * graphinter

flush(stderr()); flush(stdout())

### Name: graphinter
### Title: Graphical display of the interaction between two qualitative
###   variables
### Aliases: graphinter
### Keywords: models

### ** Examples

## Not run: 
##D data(chocolates)
##D graphinter(sensochoc, col.p = 4, col.j = 2, firstvar = 5, lastvar = 12,
##D     numr = 1, numc = 1)
## End(Not run)



cleanEx()
nameEx("hedo.cocktail")
### * hedo.cocktail

flush(stderr()); flush(stdout())

### Name: hedo.cocktail
### Title: Cocktails hedonic scores
### Aliases: hedo.cocktail
### Keywords: datasets

### ** Examples

data(cocktail)



cleanEx()
nameEx("hedochoc")
### * hedochoc

flush(stderr()); flush(stdout())

### Name: hedochoc
### Title: Chocolates hedonic scores
### Aliases: hedochoc
### Keywords: datasets

### ** Examples

data(chocolates)



cleanEx()
nameEx("histprod")
### * histprod

flush(stderr()); flush(stdout())

### Name: histprod
### Title: Histogram for each descriptor
### Aliases: histprod
### Keywords: univar

### ** Examples

data(chocolates)
histprod(sensochoc, firstvar = 5, lastvar = 10)



cleanEx()
nameEx("hsortplot")
### * hsortplot

flush(stderr()); flush(stdout())

### Name: hsortplot
### Title: Plot consumers' hierarchical sorting
### Aliases: hsortplot
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(cards)
##D group.cards<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D hsortplot(cards,group.cards)
## End(Not run)



cleanEx()
nameEx("indscal")
### * indscal

flush(stderr()); flush(stdout())

### Name: indscal
### Title: Construct the Indscal model for Napping data type
### Aliases: indscal
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
##D resindscal<- indscal(napping.don, napping.words)
##D x11()
##D prefpls(cbind(resindscal$points, napping.words))
##D x11()
##D pmfa(napping.don, napping.words, mean.conf = resindscal$points)
## End(Not run)



cleanEx()
nameEx("interact")
### * interact

flush(stderr()); flush(stdout())

### Name: interact
### Title: Estimation of interaction coefficients
### Aliases: interact
### Keywords: models

### ** Examples

## Not run: 
##D data(chocolates)
##D resinteract=interact(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
## End(Not run)



cleanEx()
nameEx("magicsort")
### * magicsort

flush(stderr()); flush(stdout())

### Name: magicsort
### Title: Returns a sorted data matrix
### Aliases: magicsort
### Keywords: manip

### ** Examples

## Example 1
data(chocolates)
resdecat<-decat(sensochoc, formul = "~Product", firstvar = 5, 
    graph = FALSE)
coltable(magicsort(resdecat$tabT), level.lower = -1.96, 
    level.upper = 1.96, main.title = "Products' description")

## Example 2
data(chocolates)
resperf<-paneliperf(sensochoc, 
    formul = "~Product+Panelist+Product:Panelist",
    formul.j = "~Product", col.j = 1, firstvar = 5, lastvar = 12,
    synthesis = FALSE, graph = FALSE)
res.sort=magicsort(resperf$prob.ind, method = "median")
coltable(res.sort, main.title = "P-values of the F-test by panelist")



cleanEx()
nameEx("nappeplot")
### * nappeplot

flush(stderr()); flush(stdout())

### Name: nappeplot
### Title: Plot panelists' tableclothe
### Aliases: nappeplot
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
## End(Not run)



cleanEx()
nameEx("nappesortplot")
### * nappesortplot

flush(stderr()); flush(stdout())

### Name: nappesortplot
### Title: Plot consumers' sorted tablecloth
### Aliases: nappesortplot
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(smoothies)
##D nappesortplot(smoothies,first="nappe")
## End(Not run)



cleanEx()
nameEx("napping")
### * napping

flush(stderr()); flush(stdout())

### Name: napping
### Title: Napping data
### Aliases: napping
### Keywords: datasets

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
##D x11()
##D pmfa(napping.don, napping.words)
## End(Not run)



cleanEx()
nameEx("napping.don")
### * napping.don

flush(stderr()); flush(stdout())

### Name: napping.don
### Title: An example of Napping data
### Aliases: napping.don
### Keywords: datasets

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
##D res <- pmfa(napping.don, napping.words)
##D res2 <- boot(napping.don,method="napping")
## End(Not run)



cleanEx()
nameEx("napping.words")
### * napping.words

flush(stderr()); flush(stdout())

### Name: napping.words
### Title: An example of "illustrative" variables to enhance results from
###   Napping data
### Aliases: napping.words
### Keywords: datasets

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
##D x11()
##D pmfa(napping.don, napping.words)
## End(Not run)



cleanEx()
nameEx("optimaldesign")
### * optimaldesign

flush(stderr()); flush(stdout())

### Name: optimaldesign
### Title: Construction of an optimal design
### Aliases: optimaldesign
### Keywords: math

### ** Examples

## Not run: 
##D optimaldesign(nbPanelist=10,nbPanelistMin=8,nbProd=5,nbProdByPanelist=3)
## End(Not run)



cleanEx()
nameEx("paneliperf")
### * paneliperf

flush(stderr()); flush(stdout())

### Name: paneliperf
### Title: Panelists' performance according to their capabilities to
###   dicriminate between products
### Aliases: paneliperf
### Keywords: models

### ** Examples

## Not run: 
##D data(chocolates)
##D res<-paneliperf(sensochoc, formul = "~Product+Panelist+Session+
##D   Product:Panelist+Product:Session+Panelist:Session",
##D   formul.j = "~Product", col.j = 1, firstvar = 5, synthesis = TRUE)
##D resprob<-magicsort(res$prob.ind, method = "median")
##D coltable(resprob, level.lower = 0.05, level.upper = 1,
##D     main.title = "P-value of the F-test (by panelist)")
##D hist(resprob,main="Histogram of the P-values",xlab="P-values")
##D 
##D resr2<-magicsort(res$r2.ind, method = "median", ascending = FALSE)
##D coltable(resr2, level.lower = 0.00, level.upper = 0.85,
##D     main.title = "Adjusted R-square (by panelist)")
##D 
##D resagree<-magicsort(res$agree, sort.mat = res$r2.ind, method = "median")
##D coltable(resagree, level.lower = 0.00, level.upper = 0.85,
##D     main.title = "Agreement between panelists")
##D hist(resagree,main="Histogram of the agreement between panelist and panel",
##D     xlab="Correlation coefficient between the product effect for 
##D     panelist and panel")
##D 
##D coltable(magicsort(res$p.value, sort.mat = res$p.value[,1], bycol = FALSE,
##D     method = "median"),
##D     main.title = "Panel performance (sorted by product P-value)")
## End(Not run)



cleanEx()
nameEx("panellipse")
### * panellipse

flush(stderr()); flush(stdout())

### Name: panellipse
### Title: Confidence ellipses around products based on panelists
###   descriptions
### Aliases: panellipse
### Keywords: multivariate

### ** Examples

## Not run: 
##D ## Example 1: PCA
##D data(chocolates)
##D res <- panellipse(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
##D coltable(res$hotelling, main.title = "P-values for the Hotelling's T2 tests")
##D 
##D ## If we consider only 12 panelists in a virtual panel, 
##D ## what would be the size of the ellipses
##D res2 <- panellipse(sensochoc, col.p = 4, col.j = 1, nbchoix = 12, firstvar = 5)
##D coltable(res2$hotelling, main.title = "P-values for the Hotelling's T2 tests")
##D 
##D ## If we want the confidence ellipses around the individual descriptions
##D panellipse(sensochoc, col.p = 4, col.j = 1, nbchoix = 1, firstvar = 5)
##D 
##D 
##D ## Example 2: MFA
##D data(chocolates)
##D res <- panellipse(sensochoc, col.p = 4, col.j = 1, firstvar = 5, 
##D     group = c(6,8), name.group = c("G1","G2"))
##D for (i in 1:dim(res$hotelling$bygroup)[3]) coltable(res$hotelling$bygroup[,,i], 
##D     main.title = paste("P-values for the Hotelling's T2 tests (",
##D     dimnames(res$hotelling$bygroup)[3][[1]][i],")",sep=""))
## End(Not run)



cleanEx()
nameEx("panellipse.session")
### * panellipse.session

flush(stderr()); flush(stdout())

### Name: panellipse.session
### Title: Repetability of panelists descriptions studied by confidence
###   ellipses around products per session
### Aliases: panellipse.session
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(chocolates)
##D res <- panellipse.session(sensochoc, col.p = 4, col.j = 1, col.s = 2, 
##D     firstvar = 5)
##D magicsort(res$variability)
##D for (i in 1:dim(res$hotelling$bysession)[3]) coltable(res$hotelling$bysession[,,i], 
##D     main.title = paste("P-values for the Hotelling's T2 tests (",
##D     dimnames(res$hotelling$bysession)[3][[1]][i],")",sep=""))
## End(Not run)



cleanEx()
nameEx("panelmatch")
### * panelmatch

flush(stderr()); flush(stdout())

### Name: panelmatch
### Title: Confidence ellipses around products based on panel descriptions
### Aliases: panelmatch
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(chocolates)
##D Panel1=sensochoc[as.numeric(sensochoc[,1])<11,]
##D Panel2=sensochoc[as.numeric(sensochoc[,1])<21 & as.numeric(sensochoc[,1])>10,]
##D Panel3=sensochoc[as.numeric(sensochoc[,1])>20,]
##D res <- panelmatch(list(P1=Panel1,P2=Panel2,P3=Panel3), col.p = 4, col.j = 1, firstvar = 5)
## End(Not run)



cleanEx()
nameEx("panelperf")
### * panelperf

flush(stderr()); flush(stdout())

### Name: panelperf
### Title: Panel's performance according to its capabilities to dicriminate
###   between products
### Aliases: panelperf
### Keywords: models

### ** Examples

data(chocolates)
res=panelperf(sensochoc, firstvar = 5, formul = "~Product+Panelist+
    Session+Product:Panelist+Session:Product+Panelist:Session")
## Sort results by product p.values.
coltable(magicsort(res$p.value, sort.mat = res$p.value[,1], bycol = FALSE,
    method = "median"), main.title = "Panel performance (sorted by product P-value)")



cleanEx()
nameEx("perfume")
### * perfume

flush(stderr()); flush(stdout())

### Name: perfume
### Title: Perfume
### Aliases: perfume
### Keywords: datasets

### ** Examples

data(perfume)

## Example of FAST
res.fast <- fast(perfume)




cleanEx()
nameEx("perfume_fcp")
### * perfume_fcp

flush(stderr()); flush(stdout())

### Name: perfume_fcp
### Title: Perfume data obtained by free choice profiling
### Aliases: perfume_fcp
### Keywords: datasets

### ** Examples

data(perfume_fcp)
res <- fcp(perfume_fcp, group = c(12,7,7,7,6,8))



cleanEx()
nameEx("plot.fahst")
### * plot.fahst

flush(stderr()); flush(stdout())

### Name: plot.fahst
### Title: Make Factorial Approach for Hierarchical Sorting Task data
###   (FAHST) graphs
### Aliases: plot.fahst
### Keywords: dplot

### ** Examples

## Not run: 
##D data(cards)
##D group<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D res.fahst <- fahst(cards,group,graph=FALSE)
##D plot.fahst(res.fahst,choix="ind",invisible="var",habillage=2,title="Cards colored according to level 2 of subject 1")
##D plot.fahst(res.fahst,choix="level",traj=TRUE)
## End(Not run)



cleanEx()
nameEx("plot.fasnt")
### * plot.fasnt

flush(stderr()); flush(stdout())

### Name: plot.fasnt
### Title: Make Factorial Approach for Sorting Napping Task data (FASNT)
###   graphs
### Aliases: plot.fasnt
### Keywords: dplot

### ** Examples

## Not run: 
##D data(smoothies)
##D res.fasnt <- fasnt(smoothies, first="nappe",graph=FALSE)
##D plot.fasnt(res.fasnt,choix="ind",invisible="var",habillage=15,
##D   title="Objects colored according to the groups provided by consumer 5")
##D plot.fasnt(res.fasnt,choix="partial",lab.partial=FALSE)
## End(Not run)



cleanEx()
nameEx("plot.fast")
### * plot.fast

flush(stderr()); flush(stdout())

### Name: plot.fast
### Title: Make Factorial Approach for Sorting Task data (FAST) graphs
### Aliases: plot.fast
### Keywords: dplot

### ** Examples

## Not run: 
##D data(perfume)
##D res.fast <- fast(perfume,graph=FALSE)
##D plot.fast(res.fast,choix="ind",invisible="var",habillage=5)
##D plot.fast(res.fast,choix="group")
## End(Not run)



cleanEx()
nameEx("plotellipse")
### * plotellipse

flush(stderr()); flush(stdout())

### Name: plotellipse
### Title: Plot confidence ellipses
### Aliases: plotellipse
### Keywords: dplot internal

### ** Examples

## Not run: 
##D data(chocolates)
##D donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
##D axe <- construct.axes(donnee, scale.unit = TRUE)
##D simul <- simulation(axe)
##D plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))
##D #######################################
##D donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
##D axe <- construct.axes(donnee, group = c(6,8), 
##D     name.group = c("A-F","T-S"),scale.unit = TRUE)
##D simul <- simulation(axe, nbgroup = (ncol(axe$partiel)-2)/(ncol(axe$moyen)-2))
##D plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))
## End(Not run)


cleanEx()
nameEx("plotpanelist")
### * plotpanelist

flush(stderr()); flush(stdout())

### Name: plotpanelist
### Title: Plotpanelist
### Aliases: plotpanelist
### Keywords: multivariate

### ** Examples

data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
plotpanelist(axe$moyen, eig = signif(axe$eig,4))



cleanEx()
nameEx("pmfa")
### * pmfa

flush(stderr()); flush(stdout())

### Name: pmfa
### Title: Procrustean Multiple Factor Analysis (PMFA)
### Aliases: pmfa
### Keywords: multivariate

### ** Examples

## Not run: 
##D data(napping)
##D nappeplot(napping.don)
##D x11()
##D pmfa(napping.don, napping.words)
## End(Not run)



cleanEx()
nameEx("print.fahst")
### * print.fahst

flush(stderr()); flush(stdout())

### Name: print.fahst
### Title: Print Factorial Approach for Hierarchical Sorting Task data
###   (FAHST) results
### Aliases: print.fahst
### Keywords: print

### ** Examples

## Not run: 
##D data(cards)
##D group<-c(2,3,3,2,2,4,2,3,2,1,3,2,3,3,3,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3)
##D res.fast <- fahst(cards,group,graph=F)
##D print.fahst(res.fahst, file="c:/fahst.csv", sep = ";")
## End(Not run)



cleanEx()
nameEx("print.fasnt")
### * print.fasnt

flush(stderr()); flush(stdout())

### Name: print.fasnt
### Title: Print Factorial Approach for Sorting Napping Task data (FASNT)
###   results
### Aliases: print.fasnt
### Keywords: print

### ** Examples

## Not run: 
##D data(smoothies)
##D res.fasnt <- fasnt(smoothies, first="nappe",graph=F)
##D print.fasnt(res.fasnt, file="c:/fasnt.csv", sep = ";")
## End(Not run)



cleanEx()
nameEx("print.fast")
### * print.fast

flush(stderr()); flush(stdout())

### Name: print.fast
### Title: Print Factorial Approach for Sorting Task data (FAST) results
### Aliases: print.fast
### Keywords: print

### ** Examples

## Not run: 
##D data(perfume)
##D res.fast <- fast(perfume,graph=FALSE)
##D print.fast(res.fast, file="c:/essai.csv", sep = ";")
## End(Not run)



cleanEx()
nameEx("scalebypanelist")
### * scalebypanelist

flush(stderr()); flush(stdout())

### Name: scalebypanelist
### Title: Scale by panelist
### Aliases: scalebypanelist
### Keywords: manip

### ** Examples

data(chocolates)
res=scalebypanelist(sensochoc, col.p = 4, col.j = 1, firstvar = 5)
res



cleanEx()
nameEx("search.desc")
### * search.desc

flush(stderr()); flush(stdout())

### Name: search.desc
### Title: Search for discriminating descriptors
### Aliases: search.desc
### Keywords: models

### ** Examples

data(chocolates)
## In this example, all the descriptos are discriminated
interesting.desc <- search.desc(sensochoc, col.j = 1, col.p = 4, 
    firstvar = 5, level = 0.5)



cleanEx()
nameEx("senso.cocktail")
### * senso.cocktail

flush(stderr()); flush(stdout())

### Name: senso.cocktail
### Title: Sensory data for 16 cocktails
### Aliases: senso.cocktail
### Keywords: datasets

### ** Examples

data(cocktail)



cleanEx()
nameEx("sensochoc")
### * sensochoc

flush(stderr()); flush(stdout())

### Name: sensochoc
### Title: Sensory data for 6 chocolates
### Aliases: sensochoc
### Keywords: datasets

### ** Examples

data(chocolates)
decat(sensochoc, formul = "~Product+Panelist", firstvar = 5, graph = FALSE)



cleanEx()
nameEx("sensopanels")
### * sensopanels

flush(stderr()); flush(stdout())

### Name: sensopanels
### Title: Sensory profiles given by 7 panels
### Aliases: sensopanels
### Keywords: datasets

### ** Examples

data(chocolates)



cleanEx()
nameEx("simulation")
### * simulation

flush(stderr()); flush(stdout())

### Name: simulation
### Title: Simulate virtual panels
### Aliases: simulation
### Keywords: models internal

### ** Examples

data(chocolates)
donnee <- cbind.data.frame(sensochoc[,c(1,4,5:18)])
axe <- construct.axes(donnee, scale.unit = TRUE)
simul <- simulation(axe)
plotellipse (simul, alpha = 0.05, eig = signif(axe$eig,4))



cleanEx()
nameEx("smoothies")
### * smoothies

flush(stderr()); flush(stdout())

### Name: smoothies
### Title: Smoothies
### Aliases: smoothies
### Keywords: datasets

### ** Examples

## Not run: 
##D data(smoothies)
##D ## Example of FASNT
##D res.fasnt<-fasnt(smoothies,first="nappe")
## End(Not run)



cleanEx()
nameEx("triangle.design")
### * triangle.design

flush(stderr()); flush(stdout())

### Name: triangle.design
### Title: Construct a design for triangle tests
### Aliases: triangle.design
### Keywords: models

### ** Examples

##Example 1
design1 = triangle.design (nbprod = 4, nbpanelist = 8)

##Example 2
design2 = triangle.design(nbprod = 4, nbpanelist = 6, bypanelist = 3,
  labprod=c("prod1","prod2","prod3","prod4"),
  labpanelist=c("John","Audrey","Peter","Martina","James","Lisa"))
  



cleanEx()
nameEx("triangle.pair.test")
### * triangle.pair.test

flush(stderr()); flush(stdout())

### Name: triangle.pair.test
### Title: Make a Triangle test for two products
### Aliases: triangle.pair.test
### Keywords: models

### ** Examples

triangle.pair.test (11, 20)  



cleanEx()
nameEx("triangle.test")
### * triangle.test

flush(stderr()); flush(stdout())

### Name: triangle.test
### Title: Make a Triangle test for a set of products
### Aliases: triangle.test
### Keywords: models

### ** Examples

design = triangle.design(nbprod = 4, nbpanelist = 6, bypanelist = 3)
answer = c("X","Y","Y","X","Z","X","Y","X","Z",
    "X","X","Z","X","Y","X","Z","X","Y")
triangle.test (design, answer)  



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
