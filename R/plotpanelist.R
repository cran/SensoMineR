"plotpanelist" <- function(mat,coord=c(1,2),name=FALSE)   {
mat=cbind.data.frame(mat[,coord],mat[,(ncol(mat)-1):ncol(mat)])
nbprod <- length(levels(mat[,ncol(mat)-1]))
lab <- mat[1:nbprod,ncol(mat)-1]
nbjuge <-  length(levels(mat[,ncol(mat)]))-1

nbpoint=nbprod*(nbjuge+1)

minx <- min(mat[1:nbpoint,1],na.rm=TRUE)
maxx <- max(mat[1:nbpoint,1],na.rm=TRUE)
miny <- min(mat[1:nbpoint,2],na.rm=TRUE)
maxy <- max(mat[1:nbpoint,2],na.rm=TRUE)

  get(getOption("device"))()
  par(mar = c(0,0,2,0))
  senso.label( mat[1:nbpoint,], clabel=0, cpoint=0,xlim=c(minx*1.05,maxx*1.05),ylim=c(1.05*miny,1.05*maxy))
  points(mat[1:nbprod,1],mat[1:nbprod,2],cex=1.2,col=1:nbprod,pch=15)
  text( mat[1:nbprod,1], mat[1:nbprod,2], labels=lab, cex = 0.8, pos = 4, offset = 0.2,col=1:nbprod)
  if (name==FALSE) for (i in 1:nbjuge) points(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],cex=0.8,col=1:nbprod,pch=20)
  if (name==TRUE) for (i in 1:nbjuge) text(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],labels=mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),ncol(mat)],cex=0.6,col=1:nbprod,pch=20)
  title(main = paste("Individual description (comp ",coord[1]," - comp ",coord[2],")",sep=""))

}
