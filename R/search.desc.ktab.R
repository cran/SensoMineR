"search.desc.ktab" <- function(ktab.donnee,level=0.5) {
  mat.desc.interessant <- ktab.donnee[[1]]
  dim.par.gr <- matrix(0,length(tab.names(ktab.donnee))-1,1)
  for (i in 1:(length(tab.names(ktab.donnee))-1) ){
  aux <- search.desc(cbind(ktab.donnee[[1]],ktab.donnee[[i+1]]),col.j=1,col.p=2,firstvar=3,level=level)
  dim.par.gr[i]<- ncol(aux)-2               # on enleve la colonne juge et produit
  if (dim.par.gr[i]==0) print(paste("There's no significant descriptors in the group of variables",i))
  mat.desc.interessant <- cbind(mat.desc.interessant,aux[,-c(1,2)])
  }
  res.ktab <- ktab.data.frame(mat.desc.interessant,blocks=c(2,dim.par.gr),tabnames=tab.names(ktab.donnee))
  return(res.ktab)
}
