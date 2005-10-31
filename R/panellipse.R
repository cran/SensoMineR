"panellipse" <- function(donnee,col.p,col.j,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,bloc=NULL,name.bloc=NULL,level.search.desc=0.5,centerbypanelist=TRUE,scalebypanelist=FALSE,name.panelist=FALSE,cex=1,color=NULL){

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

  if (length(color)==0) color = c("black","red","green3","blue",
    "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
    "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
    "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
    "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")

  don <- cbind.data.frame(donnee[,col.j],donnee[,col.p],donnee[,firstvar:lastvar])
  colnames(don) <- colnames(donnee)[c(col.j,col.p,firstvar:lastvar)]
  if (length(bloc)<2) ktab.donnee <- ktab.data.frame(don,blocks=c(2,lastvar-firstvar+1),tabnames=c("JP","Gr1"))
  if (length(bloc)>1) ktab.donnee <- ktab.data.frame(don,blocks=c(2,bloc),tabnames=c("JP",name.bloc))
  ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=level.search.desc)
  axe <- construct.axes(ktab.interesting.desc,coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist)
  labprod=axe$moyen[axe$moyen[,ncol(axe$moyen)]==0,ncol(axe$moyen)-1]

  nbprod = length(labprod)

if (length(bloc)<2) {
  mat = matrix(NA,nbprod,nbprod)
    aa = axe$moyen[-(1:nbprod),]
    for (i in 1:(nbprod-1)){
      for (j in i:nbprod){
        if (length(nbchoix)==0) mat[i,j] = mat[j,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord])
        if (length(nbchoix)!=0) mat[i,j] = mat[j,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
    }}
    diag(mat)=1
    colnames(mat)=rownames(mat)=labprod
}

if (length(bloc)>1) {
  mat = array(NA,dim=c(nbprod,nbprod,length(bloc)+1))
  for (k in 1:length(bloc)){
    aa = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(k-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
    for (i in 1:(nbprod-1)){
      for (j in i:nbprod){
        if (length(nbchoix)==0) mat[i,j,k] = mat[j,i,k] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord])
        if (length(nbchoix)!=0) mat[i,j,k] = mat[j,i,k] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
    }}
    for (i in 1:nbprod) mat[i,i,k]=1
  }
  aa=axe$moyen[-(1:length(labprod)),]
  for (i in 1:(length(labprod)-1)){
    for (j in (i+1):length(labprod)){
      if (length(nbchoix)==0) mat[i,j,length(bloc)+1] = mat[j,i,length(bloc)+1] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord])
      if (length(nbchoix)!=0) mat[i,j,length(bloc)+1] = mat[j,i,length(bloc)+1] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
  }}
  for (i in 1:nbprod) mat[i,i,length(bloc)+1]=1
  dimnames(mat)=list(labprod,labprod,c(paste("Group",1:length(bloc),sep=" "),"global"))

  mat2 = array(NA,dim=c(length(bloc),length(bloc),nbprod))
  for (i in 1:nbprod){
    for (k in 1:(length(bloc)-1)){
      aa = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(k-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
      for (kk in (k+1):length(bloc)){
        aa2 = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(kk-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
        if (length(nbchoix)==0) mat2[k,kk,i] = mat2[kk,k,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa2[aa2[,ncol(aa2)-1]==labprod[i],coord])
        if (length(nbchoix)!=0) mat2[k,kk,i] = mat2[kk,k,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa2[aa2[,ncol(aa2)-1]==labprod[i],coord],nbchoix,nbchoix)
      }}
    for (k in 1:length(bloc)) mat2[k,k,i]=1
  }
  dimnames(mat2)=list(c(paste("Group",1:length(bloc),sep=" ")),c(paste("Group",1:length(bloc),sep=" ")),labprod)
}
  
  plotpanelist(axe$moyen,coord=coord,eig=axe$eig,color=color,name=name.panelist,cex=cex)
  get(getOption("device"))(12,8)
  simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbchoix=nbchoix,nbsimul=nbsimul)
  plotellipse(simul,alpha=alpha,coord=coord,eig=axe$eig,color=color,cex=cex)  
  res <- list()
  res$eig= axe[[length(names(axe))]]
  res$coordinates= axe[-length(names(axe))]
  if (length(bloc)<2) res$hotelling=mat
  if (length(bloc)>1) res$hotelling=list(bygroup=mat,byproduct=mat2)
  return(res)
}
