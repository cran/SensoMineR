"plotellipseinterhmfa" <- function(mat,alpha=0.05,coord=c(1,2),nbbloc=1,eig,cex=1,color=NULL,hmfa){

#################################################################
"ellipse" <- function(loc, cov,alpha)
      {
            A <- cov
           detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
            dist <- sqrt(qchisq(1-alpha/2, 2))
            ylimit <- sqrt(A[2, 2]) * dist
            y <- seq( - ylimit, ylimit, 0.01 * ylimit)
            sqrt.discr <- sqrt(detA/A[2, 2]^2 * abs(A[2, 2] * dist^2 - y^2))
            sqrt.discr[c(1, length(sqrt.discr))] <- 0
            b <- loc[1] + A[1, 2]/A[2, 2] * y
            x1 <- b - sqrt.discr
            x2 <- b + sqrt.discr
            y <- loc[2] + y
            return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
      }
#################################################################

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")

  matmoyP=cbind.data.frame(mat$moy$P[,coord],mat$moy$P[,ncol(mat$moy$P)])
  matmoyPJ=cbind.data.frame(mat$moy$PJ[,coord],mat$moy$PJ[,ncol(mat$moy$PJ)])
  matmoysimul=cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
  matP=cbind.data.frame(mat$partiel$P[,coord],mat$partiel$P[,ncol(mat$partiel$P)])
  matPJ=cbind.data.frame(mat$partiel$PJ[,coord],mat$partiel$PJ[,ncol(mat$partiel$PJ)])
  matsimul=cbind.data.frame(mat$partiel$simul[,coord],mat$partiel$simul[,ncol(mat$partiel$simul)])

hierar=hmfa[[1]]
ptpartiel = hmfa[[2]]
nbp <- nrow(matP)
nbprod <- nbp/nbbloc
coord.ellipse.a.tracer <- matrix(0,402,2*nbp)

p <- 2
nbjuge <-  nrow(matPJ)/nrow(matP)
nbsimul <-  nrow(matsimul)/nrow(matP)

for (i in 1:nbp){
  VX <- var(matsimul[((i-1)*nbsimul+1):(i*nbsimul),1:2])
  coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse(t(matP[i,1:2]),VX,alpha)
}

minx <- min(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
maxx <- max(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
miny <- min(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)
maxy <- max(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)
  par(mar = c(0,0,2,0))
  senso.label( matP[,1:2], clabel=0, cpoint=0,xlim=c(minx*1.05,maxx*1.05),ylim=c(1.05*miny,1.05*maxy))
  points(matmoyP[,1], matmoyP[,2],cex=0.9*cex,col=color[1:nbprod],pch=15)
  for (j in 1:length(ptpartiel[[length(hierar)]])){
    for (i in 1:nbprod) {
      points(ptpartiel[[length(hierar)]][[j]][,1], ptpartiel[[length(hierar)]][[j]][,2],cex=0.8*cex,col=color[1:nbprod],pch=20)
      lines(c(ptpartiel[[length(hierar)]][[j]][i,1],matmoyP[i,1]),c( ptpartiel[[length(hierar)]][[j]][i,2],matmoyP[i,2]),col=color[i],lwd=1.5)
  }}
  text( matmoyP[,1], matmoyP[,2], matmoyP[,ncol(matmoyP)], cex = 0.8*cex, pos = 4, offset = 0.2,col=color[1:nbprod])
  for (j in 1:nbbloc){
    for (i in 1:nbprod) {
      points(ptpartiel[[length(hierar)-1]][[j]][i,1], ptpartiel[[length(hierar)-1]][[j]][i,2],cex=0.8*cex,col=color[i],pch=20)
      jj=1
      for (k in 2: length(hierar[[length(hierar)]])){
        if (j > sum(hierar[[length(hierar)]][1:(k-1)])) jj=jj+1
      }
      lines(c(ptpartiel[[length(hierar)-1]][[j]][i,1],ptpartiel[[length(hierar)]][[jj]][i,1]),c( ptpartiel[[length(hierar)-1]][[j]][i,2],ptpartiel[[length(hierar)]][[jj]][i,2]),col=color[i],lty=j)
      lines(coord.ellipse.a.tracer[,(1+2*((i+(j-1)*nbprod)-1)):(2*(i+(j-1)*nbprod))],col=color[i],lty=j)
    }
  }
  title(main = paste("Confidence ellipses for the partial points: comp ",coord[1]," (",eig[coord[1],2],"%) - comp ",coord[2]," (",eig[coord[2],2],"%)",sep=""))
}
