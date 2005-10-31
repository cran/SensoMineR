"graphinter" <- function(donnee, col.p, col.j, firstvar, lastvar=ncol(donnee), numr=2,numc=2) {

     for (j in 1:(firstvar-1)) donnee[,j] <- as.factor(donnee[,j])
     nbprod <- length(levels(donnee[,col.p]))
     nbdesc <- lastvar-firstvar+1
     nbseance <- length(levels(donnee[,col.j]))
          
     moy <- matrix(0,nbprod,nbdesc)
     for (p in 1:nbprod) moy[p,] <- mean(donnee[(donnee[,col.p]==levels(donnee[,col.p])[p]),(firstvar:lastvar)])
     dimnames(moy) <- list(levels(donnee[,col.p]),variable.names(donnee)[firstvar:lastvar])

     moyS <- array(0,c(nbprod,nbdesc,nbseance))
     for (p in 1: nbprod) {
         for (s in 1:nbseance) {
           moyS[p,,s] <- mean(donnee[((donnee[,col.p]==levels(donnee[,col.p])[p])&(donnee[,col.j]==levels(donnee[,col.j])[s])),(firstvar:lastvar)]) 
         }
     }
     dimnames(moyS) <- list(levels(donnee[,col.p]),variable.names(donnee)[firstvar:lastvar],levels(donnee[,col.j]))

mult <- nbdesc %/% (numr*numc)
if (nbdesc==(nbdesc %/% (numr*numc))*(numr*numc)) mult=mult-1
for (m in 0:mult) {
    par(mfrow = c(numr,numc))

    for (nbd in 1:(numr*numc)) {
          nb <- (m*(numr*numc)+nbd)
          if (nb <= nbdesc) {
              xmin = ymin <- min(moy[,nb],moyS[,nb,])-0.2
              xmax = ymax <- max(moy[,nb],moyS[,nb,])+0.2
        for (s in 1:nbseance) {
                if (s==1) {
            plot(sort(moy[,nb]),sort(moyS[,nb,s]),type="o",xlab=paste("Mean on the whole ",colnames(donnee)[col.j],"s",sep=""), ylab=paste("Mean per ",colnames(donnee)[col.j],sep=""),
                     cex.lab = 0.8, asp = 1, pch = 20, xlim = c(xmin,xmax), ylim = c(ymin,ymax), col = "violetred4")
                    senso.scatterutil.eti(sort(moy[,nb]),sort(moyS[,nb,s]),label=levels(donnee[,col.p]),clabel=1, pos = 3,offset = 0.4, font = 1)
                }
            else {
            points(sort(moy[,nb]),sort(moyS[,nb,s]),type="o",pch=20,col=s)
                    senso.scatterutil.eti(sort(moy[,nb]),sort(moyS[,nb,s]),label="",clabel=0,coul=rep(s,length(moy)))
                title(variable.names(moy)[nb])
            }
        }
          } 
    }
if (m < mult) get(getOption("device"))()
}
moyenne <- list()
moyenne$col.p <- moy
moyenne$col.j <- moyS
return(moyenne)
}
