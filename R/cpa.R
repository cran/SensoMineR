cpa<- function(senso, hedo, coord=c(1,2),center=TRUE,scale=TRUE,nb.clusters=0,scale.unit=FALSE,col=terrain.colors(45)[1:41]) {

    if (max(coord) > (nrow(hedo)-1)) {
      print (paste("Problem with coord. Max (coord) must be less than",nrow(hedo)-1," Axes 1-2 will be taken",sep=""))
      coord=c(1,2)
    }
    senso <- scale(senso,center=center,scale=scale)[,]
    hedo <- scale(hedo,center=center,scale=scale)[,]
    if (scale) senso=senso*sqrt(nrow(senso)/(nrow(senso)-1))
    if (scale) hedo=hedo*sqrt(nrow(hedo)/(nrow(hedo)-1))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    senso <- data.frame(senso)
    hedo <- data.frame(hedo)
    nbjuge <- ncol(hedo)
    nbdesc <- ncol(senso)

  hc <- hclust(dist(t(hedo)),method="ward")
  plot(as.dendrogram(hc),main="Cluster Dendrogram",xlab="Panelists",leaflab="none") 
  get(getOption("device"))()
  if (nb.clusters==0){
    classif=hopach(t(hedo),d="euclid",K=10,mss="mean")
    nb.clusters=classif$clustering$k
  }
  clusters=pam(t(hedo),k=nb.clusters)$clustering
  mat <- matrix(0,nb.clusters,nrow(hedo))
  dimnames(mat) <- list(1:nb.clusters,rownames(hedo))
  for (i in 1:nb.clusters){
    mat[i,] <- apply(t(hedo[,clusters==i]),2,mean)
    rownames(mat)[i] <- paste("cluster",i)
  }  
  desc.clusters=cor(senso,t(mat))

    A <- rbind.data.frame(t(hedo),mat,t(senso))
    colnames(A) <- row.names(hedo)
    result <- A

    hedo.pca <- pca(A,supind=(nbjuge+1):nrow(A),coord=coord,scale.unit=scale.unit,clabel=0.7,graph=TRUE)
    get(getOption("device"))()
    par(mar = c(0,0,2,0))
    par(mar = c(0,0,2,0))
    senso.label(hedo.pca$li[,coord], clabel=0, cpoint=0)
    points(hedo.pca$li[,coord], pch=20,cex=0.8,col=2+clusters)
    text(hedo.pca$li[,coord[1]],hedo.pca$li[,coord[2]], labels = row.names(A[1:nbjuge,]),col=2+clusters, cex = 0.7, pos = 4, offset = 0.2)
    title(main = paste("Individuals factor map: comp",coord[1]," (",signif(100*hedo.pca$eig[coord[1]]/sum(hedo.pca$eig),4),"%)","-comp",coord[2]," (",signif(100*hedo.pca$eig[coord[2]]/sum(hedo.pca$eig),4),"%)",sep=""),cex.main = 1.1, font.main = 2, col.main = 1, adj = 0.5)

    points(hedo.pca$lisup[1:nb.clusters,coord], col=(3:(2+nb.clusters)),pch=15,cex=0.9)
    text(hedo.pca$lisup[1:nb.clusters,coord[1]],hedo.pca$lisup[1:nb.clusters,coord[2]], labels = row.names(A[(nbjuge+1):(nbjuge+nb.clusters),]), cex = 0.9, pos = 1, offset = 0.05, col = (3:(2+nb.clusters)))
    points(hedo.pca$lisup[(nb.clusters+1):nrow(hedo.pca$lisup),coord], col="red",pch=20,cex=0.8)
    text(hedo.pca$lisup[(nb.clusters+1):nrow(hedo.pca$lisup),coord[1]],hedo.pca$lisup[(nb.clusters+1):nrow(hedo.pca$lisup),coord[2]], labels = row.names(A[(nbjuge+nb.clusters+1):nrow(A),]), cex = 0.7, pos = 1, offset = 0.05, col = "red")

#    senso.sup.classif <- suprow(hedo.pca,mat)
#    senso.sup <- suprow(hedo.pca,t(senso))
#    points(senso.sup.classif$lisup[,coord], col=(3:(2+nb.clusters)),pch=15,cex=0.9)
#    text(senso.sup.classif$lisup[,coord[1]],senso.sup.classif$lisup[,coord[2]], labels = row.names(A[(nbjuge+1):(nbjuge+nb.clusters),]), cex = 0.9, pos = 1, offset = 0.05, col = (3:(2+nb.clusters)))
#    points(senso.sup$lisup[,coord], col="red",pch=20,cex=0.8)
#    text(senso.sup$lisup[,coord[1]],senso.sup$lisup[,coord[2]], labels = row.names(A[(nbjuge+nb.clusters+1):nrow(A),]), cex = 0.7, pos = 1, offset = 0.05, col = "red")

    TA <- t(A)
    coef <- matrix(0,nbjuge+nb.clusters,nbdesc)
    for (d in 1:nbdesc) {
      coef[1:nbjuge,d] <- cor(TA[,1:nbjuge],TA[,nbjuge+nb.clusters+d])
      coef[(nbjuge+1):(nbjuge+nb.clusters),d] <- cor(TA[,(nbjuge+1):(nbjuge+nb.clusters)],TA[,nbjuge+nb.clusters+d])
    }
    coef <- data.frame(coef)
    colnames(coef) <- colnames(senso)
    B <- cbind.data.frame(rbind(hedo.pca$li,hedo.pca$lisup[1:nb.clusters,]),coef)
    for (d in 1:nbdesc) {
      get(getOption("device"))(9,7)
      par(mar = c(4.2,4.1,3.5,2))
      colplot(B, k=nb.clusters,coord, (nrow(hedo)+d),col=col, xlab=paste("Comp",coord[1]," (",signif(100*hedo.pca$eig[coord[1]]/sum(hedo.pca$eig),4),"%)",sep=""), ylab=paste("Comp",coord[2]," (",signif(100*hedo.pca$eig[coord[2]]/sum(hedo.pca$eig),4),"%)",sep=""))
      points(hedo.pca$lisup[nb.clusters+d,coord[1]],hedo.pca$lisup[nb.clusters+d,coord[2]],col="red",pch=15,cex=0.8)
      text(hedo.pca$lisup[nb.clusters+d,coord[1]],hedo.pca$lisup[nb.clusters+d,coord[2]],col="red",labels=colnames(B)[nrow(hedo)+d],pos = 1, offset = 0.05)
      title(main = paste("Consumers' preferences analysed by",colnames(B)[nrow(hedo)+d]),cex.main = 1.1, font.main = 2)
    }
    don <- cbind.data.frame(as.factor(clusters),t(hedo))
    colnames(don) <- c("clusters",rownames(hedo))
    resdecat <- decat(don,formul="~clusters",firstvar=2,proba=1,graph=FALSE)
    res <- list()
    res$clusters <- clusters
    res$result <- result
    res$prod.clusters <- resdecat$resT
    res$desc.clusters <- desc.clusters
    return(res)
}

colplot<-function(mat, k=0,coord, z, level=41, col = terrain.colors(level+level%/%10)[1:level], xlab="", ylab="") { #heat.colors(level)

    abs <- coord[1]
    ord <- coord[2]
    x <- mat[,abs]
    y <- mat[,ord]
    z <- mat[,z]

#    x1 <- min(z)
#    x2 <- max(z)
    x1 <- -1
    x2 <- 1

    plot(mat[,abs],mat[,ord],xlab=xlab, ylab=ylab,asp=1,type="n")
    legend(min(x)*1.15,max(y),c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1),fill=c(col[level],col[(level%/%2)+(level%/%4)+(level%/%8)+1],col[(level%/%2)+(level%/%4)+1],col[(level%/%2)+(level%/%8)+1],col[(level%/%2)+1],col[(level%/%4)+(level%/%8)+1],col[(level%/%4)+1],col[(level%/%8)+1],col[1]),cex=0.7)
    abline(v=0,lty=2)
    abline(h=0,lty=2)
#    senso.label(mat[,c(abs,ord)],cpoint=0,clabel=0)
#    legend(min(x),max(y),c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1),fill=c(col[level],col[(level%/%2)+(level%/%4)+(level%/%8)+1],col[(level%/%2)+(level%/%4)+1],col[(level%/%2)+(level%/%8)+1],col[(level%/%2)+1],col[(level%/%4)+(level%/%8)+1],col[(level%/%4)+1],col[(level%/%8)+1],col[1]),cex=0.7)
####rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    n <- nrow(mat)
    h <- (x2-x1)/level

    for (ind in 1:(n-k)) points(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],pch=20)
    for (ind in (n-k+1):n) points(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],pch=15,cex=1)
    for (ind in (n-k+1):n) text(x[ind],y[ind],col=col[max(1,(z[ind]-x1)%/%h)],rownames(mat)[ind],cex=1,pos = 1, offset = 0.05)
}
