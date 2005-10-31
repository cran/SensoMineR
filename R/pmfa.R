pmfa<-function(matrice,matrice.illu=NULL,mean.conf=NULL,dilat=TRUE,graph.ind=TRUE,lim=c(60,40),coord=c(1,2))  {

coeffRV<-function(X,Y){
X <- scale(X,scale=FALSE)
Y <- scale(Y,scale=FALSE)
W1<-X%*%t(X)
W2<-Y%*%t(Y)
rv<-sum(diag(W1%*%W2))/(sum(diag(W1%*%W1))*sum(diag(W2%*%W2)))^0.5
return(rv)
}

procrustes <- function(amat, target, orthogonal =FALSE, translate =FALSE, magnify =FALSE)
{
for (i in nrow(amat):1){   # à l'envers pour pas gérer les lignes effacées
  if (any(is.na(amat)[i,])|any(is.na(target)[i,])) {
    amat <- amat [-i,]
    target <- target[-i,]
  }
}  
      dA <- dim(amat)
      dX <- dim(target)
      if(length(dA) != 2 || length(dX) != 2)
            stop("arguments amat and target must be matrices")
      if(any(dA != dX))
            stop("dimensions of amat and target must match")
      if(length(attr(amat, "tmat")))
            stop("oblique loadings matrix not allowed for amat")
      if(orthogonal) {
            if(translate) {
                  p <- dX[1]
                  target.m <- (rep(1/p, p) %*% target)[,  ]
                  amat.m <- (rep(1/p, p) %*% amat)[,  ]
                  target.c <- scale(target, center = target.m, scale =FALSE)
                  amat.c <- scale(amat, center = amat.m, scale =FALSE)
                  j <- svd(crossprod(target.c, amat.c))
            }
            else {
                  amat.c <- amat
                  j <- svd(crossprod(target, amat))
            }
            rot <- j$v %*% t(j$u)
            if(magnify)
                  beta <- sum(j$d)/sum(amat.c^2)
            else beta <- 1
            B <- beta * amat.c %*% rot
            if(translate)
                  B <- B + rep(as.vector(target.m), rep.int(p, dX[2]))
            value <- list(rmat = B, tmat = rot, magnify = beta)
            if(translate)
                  value$translate <- target.m - (rot %*% amat.m)[,  ]
      }
      else {
            b <- solve(amat, target)
            gamma <- sqrt(diag(solve(crossprod(b))))
            rot <- b * rep(gamma, rep.int(dim(b)[1], length(gamma)))
            B <- amat %*% rot
            fcor <- solve(crossprod(rot))
            value <- list(rmat = B, tmat = rot, correlation = fcor)
      }
      value
}

nbjuge <- ncol(matrice)/2
matrice <- scale(matrice,center=TRUE,scale=FALSE)
do.mfa=FALSE
if (length(mean.conf)==0){
  do.mfa=TRUE
  res.afm <- afmult(ktab.data.frame(as.data.frame(matrice),blocks=rep(2,nbjuge)),scann=FALSE,nf=max(coord))
  mean.conf <- res.afm$li
}
if (!do.mfa){
  aa=cor(matrice,mean.conf[,coord],use="pairwise.complete.obs")
  senso.corcircle(aa, fullcircle=TRUE)
  title(main = paste("Correlation circle (comp",coord[1],"-","comp",coord[2],")"))
  get(getOption("device"))()

  par(mar = c(0,0,2,0))
  senso.label( mean.conf[,coord], clabel=0, cpoint=0.8, include.origin = FALSE)
  text(mean.conf[,coord[1]], mean.conf[,coord[2]], labels = rownames(mean.conf), cex = 0.8, pos = 4, offset = 0.2)
  title(main = paste("Individuals factor map (","comp",coord[1],"-","comp",coord[2],")"))
  get(getOption("device"))()
}

if (length(matrice.illu)!=0){
  bb=cor(matrice.illu,mean.conf[,coord])
  senso.corcircle(bb, fullcircle=TRUE,col="blue")
  title(main = paste("Correlation circle (comp",coord[1],"-","comp",coord[2],")"))
  get(getOption("device"))()
}

mean.conf <- as.matrix(mean.conf[,coord])
res <- matrix(0,nbjuge,1)
for (j in 1:nbjuge){
  atourner <- as.matrix(matrice[,(2*(j-1)+1):(2*j)])
  if ((dilat==TRUE)&(do.mfa==TRUE)){
    eig <- eigen(1/ncol(atourner)*t(atourner)%*%atourner)
    res.procrustes<-procrustes(atourner, mean.conf, orthogonal=TRUE, translate=TRUE, magnify=FALSE)
    magnify <- res.afm$eig[1]/sqrt(eig$values[1])
  }
  else {
    res.procrustes<-procrustes(atourner, mean.conf, orthogonal=TRUE, translate=TRUE, magnify=TRUE)
    magnify <- res.procrustes$magnify
  }
  tourne <- atourner%*%res.procrustes$tmat * magnify
  res[j] <- coeffRV(mean.conf,tourne)

  if (graph.ind==TRUE){
    dd=cbind(mean.conf,tourne)
    nappe <- rbind(matrix(c(-lim[1]/2,-lim[2]/2,-lim[1]/2,lim[2]/2,lim[1]/2,lim[2]/2,lim[1]/2,-lim[2]/2),ncol=2,byrow=TRUE)%*%res.procrustes$tmat * magnify)
  
    if (j!=1) get(getOption("device"))()
    plot(rbind(tourne,mean.conf,nappe),type="n",xlab=paste("Dim",coord[1]),ylab=paste("Dim",coord[2]),asp=1,main=colnames(matrice)[2*j],sub=paste("RV between the mean representation and the representation of",colnames(matrice)[2*j],": ",signif(res[j],4)),cex.sub=0.8)
    for (i in 1:nrow(mean.conf)) points(tourne[i,1],tourne[i,2],cex=0.8,pch=20,col=3)
    for (i in 1:nrow(mean.conf)) text(tourne[i,1],tourne[i,2],rownames(matrice)[i],col=3,font=3,cex=0.8, pos = 4, offset = 0.2)
    for (i in 1:nrow(mean.conf)) points(mean.conf[i,1],mean.conf[i,2],cex=0.8,pch=20)
    for (i in 1:nrow(mean.conf)) text(mean.conf[i,1],mean.conf[i,2],rownames(matrice)[i],cex=0.8, pos = 4, offset = 0.2)
    lines(nappe[c(1:4,1),],col=3,lty=2)
  }
 }
dimnames(res) <- list(colnames(matrice)[(1:(ncol(matrice)/2))*2],"RV coeff")
return(res)
}
