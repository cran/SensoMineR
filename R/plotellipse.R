"plotellipse" <- function(mat,alpha=0.05,coord=c(1,2)){

plotellipseinter(mat,alpha=alpha,coord=coord,nbbloc=1,moy=TRUE)
if (length(mat$partiel)!=0) {
  get(getOption("device"))()
  nbbloc=length(levels(mat$partiel[,ncol(mat$partiel)])) / length(levels(mat$moy[,ncol(mat$moy)]))
  plotellipseinter(mat,alpha=alpha,coord=coord,nbbloc=nbbloc,moy=FALSE)
}
}
