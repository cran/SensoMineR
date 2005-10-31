"search.desc" <- function(matrice,col.j,col.p,firstvar,lastvar=ncol(matrice),level=0.5){
  lab<-labels(matrice)[[2]]
  nomdescripteur<-lab[firstvar:lastvar]
  for (i in 1:(firstvar-1)) matrice[,i]<-as.factor(matrice[,i])
  tab.F<-matrix(0,(lastvar-firstvar+1),1)
  for (i in firstvar:lastvar){
    aux <- summary(aov(as.formula(paste(lab[i],"~",lab[col.p],"+",lab[col.j])),data=matrice,na.action=na.exclude))
    tab.F[i-firstvar+1] <- pf(aux[[1]][1,4],aux[[1]][1,1],aux[[1]][(dim(aux[[1]])[[1]]),1],lower.tail=FALSE)
  }
  dimnames(tab.F)<-list(nomdescripteur,NULL)
  resF <- vector("list",length=1)
  select <- (1:nrow(tab.F))
  resF <- data.frame(Variables=as.factor(dimnames(tab.F)[[1]][rev(order(tab.F))][select]),Proba=as.numeric(tab.F[rev(order(tab.F))][select]))
  mat.analyse <- data.frame(as.factor(matrice[,1]))
  for (i in 2:(firstvar-1)) mat.analyse <- cbind.data.frame(mat.analyse,matrice[,i])
  for (i in firstvar:lastvar){
    if (tab.F[i-firstvar+1]<level) mat.analyse <- cbind(mat.analyse,matrice[,i])
  }
  dimnames(mat.analyse)[[2]][1:(firstvar-1)]<-lab[1:(firstvar-1)]
  dimnames(mat.analyse)[[2]][firstvar:dim(mat.analyse)[2]]<-dimnames(tab.F)[[1]][tab.F<level]
  return(mat.analyse)
}
