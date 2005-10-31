"plotellipse" <- function(mat,alpha=0.05,coord=c(1,2),eig){

plotellipseinter(mat,alpha=alpha,coord=coord,nbbloc=1,moy=TRUE,eig=eig)
if (length(mat$partiel)!=0) {
  get(getOption("device"))(max(14,eig[coord[1],1]/eig[coord[1],2]*8),8)
  nbbloc=length(levels(mat$partiel[,ncol(mat$partiel)])) / length(levels(mat$moy[,ncol(mat$moy)]))
  plotellipseinter(mat,alpha=alpha,coord=coord,nbbloc=nbbloc,moy=FALSE,eig=eig)
}
}
