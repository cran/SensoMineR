"panellipse" <- function(donnee,col.p,col.j,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,bloc=NULL,name.bloc=NULL,level.search.desc=0.5,centerbypanelist=TRUE,scalebypanelist=FALSE){

hotelling <- function(d1,d2,n1=nrow(d1),n2=nrow(d2)){
    k <- ncol(d1)
#    n1 <- nrow(d1)
#    n2 <- nrow(d2)
    xbar1 <- apply(d1,2,mean)
    xbar2 <- apply(d2,2,mean)
    dbar <- xbar2-xbar1
    v <- ((n1-1)*var(d1)+(n2-1)*var(d2))/(n1+n2-2)
    t2 <- n1*n2*dbar%*%solve(v)%*%dbar/(n1+n2)
    f <- (n1+n2-k-1)*t2/((n1+n2-2)*k)
    return(pf(f,k,n1+n2-k-1,lower.tail=FALSE))
}

  don <- cbind.data.frame(donnee[,col.j],donnee[,col.p],donnee[,firstvar:lastvar])
  colnames(don) <- colnames(donnee)[c(col.j,col.p,firstvar:lastvar)]
  if (length(bloc)<2) ktab.donnee <- ktab.data.frame(don,blocks=c(2,lastvar-firstvar+1),tabnames=c("JP","Gr1"))
  if (length(bloc)>1) ktab.donnee <- ktab.data.frame(don,blocks=c(2,bloc),tabnames=c("JP",name.bloc))
  ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=level.search.desc)
  axe <- construct.axes(ktab.interesting.desc,coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist)
  labprod=axe$moyen[axe$moyen[,ncol(axe$moyen)]==0,ncol(axe$moyen)-1]
  aa=axe$moyen[-(1:length(labprod)),]
  mat = matrix(NA,length(labprod),length(labprod))
  for (i in 1:(length(labprod)-1)){
    for (j in i:length(labprod)){
      if (length(nbchoix)==0) mat[i,j] = mat[j,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord])
      if (length(nbchoix)!=0) mat[i,j] = mat[j,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
  }}
  diag(mat)=1
  colnames(mat)=rownames(mat)=labprod
  plotpanelist(axe$moyen,coord=coord,eig=axe$eig)
  get(getOption("device"))(12,8)
  simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbchoix=nbchoix,nbsimul=nbsimul)
  plotellipse(simul,alpha=alpha,coord=coord,eig=axe$eig)  
  res <- list()
  res$eig= axe[[length(names(axe))]]
  res$coordinates= axe[-length(names(axe))]
  res$hotelling=mat
  return(res)
}
