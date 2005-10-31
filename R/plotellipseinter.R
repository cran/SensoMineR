"plotellipseinter" <- function(mat,alpha=0.05,coord=c(1,2),nbbloc=1,moy=TRUE){

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

if (moy ==TRUE) mat=cbind.data.frame(mat$moy[,coord],mat$moy[,ncol(mat$moy)])
if (moy ==FALSE) {
  matmoy <- cbind.data.frame(mat$moy[,coord],mat$moy[,ncol(mat$moy)])
  mat=cbind.data.frame(mat$partiel[,coord],mat$partiel[,ncol(mat$partiel)])
}  
nbp <- length(levels(mat[,ncol(mat)]))
nbprod <- nbp/nbbloc
X <- matrix(0,2,nbp)
lab <- mat[1:nbp,ncol(mat)]
VX <- matrix(0,2,2)
coord.ellipse.a.tracer <- matrix(0,402,2*nbp)

p <- 2
nbjuge <-  which(mat[,ncol(mat)]==lab[2])[2] -nbp -1
nbsimul <-  nrow(mat) - which(mat[,ncol(mat)]==lab[nbp-1])[length(which(mat[,ncol(mat)]==lab[nbp-1]))]
for (i in 1:nbp){
  fin <- dim(mat)[[1]]-nbsimul*(nbp-i)
  X[,i] <- t(apply(mat[(fin-nbsimul+1):fin,1:2],2,mean))
  VX <- var(mat[(fin-nbsimul+1):fin,1:2])
  coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse(X[,i],VX,alpha)
}
minx <- min(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
maxx <- max(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
miny <- min(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)
maxy <- max(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)

  par(mar = c(0,0,2,0))
  senso.label( t(X), clabel=0, cpoint=0,xlim=c(minx*1.05,maxx*1.05),ylim=c(1.05*miny,1.05*maxy))
  if (moy==FALSE) points(matmoy[1:nbprod,1], matmoy[1:nbprod,2],cex=0.8,col=1:nbprod,pch=15)
  for (j in 1:nbbloc){
    text( t(X)[,1], t(X)[,2], lab, cex = 0.8, pos = 4, offset = 0.2,col=1:nbprod)
    for (i in 1:nbprod) {
      points(t(X)[i+(j-1)*nbprod,1], t(X)[i+(j-1)*nbprod,2],cex=0.8,col=i,pch=20)
      if (moy==FALSE) lines(c(t(X)[i+(j-1)*nbprod,1],matmoy[i,1]),c( t(X)[i+(j-1)*nbprod,2],matmoy[i,2]),col=i,lty=j)
    }
    for (i in 1:nbprod) lines(coord.ellipse.a.tracer[,(1+2*((i+(j-1)*nbprod)-1)):(2*(i+(j-1)*nbprod))],col=i,lty=j)
  }
  if (moy==TRUE) title(main = paste("Confidence ellipses for the mean points (comp ",coord[1]," - comp ",coord[2],")",sep=""))
  if (moy==FALSE) title(main = paste("Confidence ellipses for the partial points (comp ",coord[1]," - comp ",coord[2],")",sep=""))
}
