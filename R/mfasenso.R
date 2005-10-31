############### Programme d'AF : 2 ktableaux le 1er actif, le 2ème illustratif
############## ressort les coord des points moyens et partiels, des individus actifs et illustratifs
######## Si dans les k-tableaux il n'y a qu'un groupe, fait aussi des ACP.

mfasenso <- function(ktab,ktab.illu=NULL,scale.unit=TRUE,nbcoord=2,poids=NULL){

##########################################################################
nb.column <- function(ktab,num.group){
  nb.col <- list()
  auxi <- 0
  auxi2 <- 0
  if (num.group==1) nb.col$prec <- 0
  else {
    for (i in 1:(num.group-1)) auxi <- auxi + ktab$blo[i]
  }
  if (num.group==length(ktab$blo)) nb.col$suiv <- 0
  else {
    for (i in (num.group+1):length(ktab$blo)) auxi2 <- auxi2 + ktab$blo[i]
  }
  nb.col$prec <- auxi
  nb.col$suiv <- auxi2
  return(nb.col)
}
##########################################################################

  if (length(poids)==0) poids=rep(1,length(tab.names(ktab)))
  tab <- matrix(0,nrow(ktab[[1]]),0)
  for (i in 1:length(tab.names(ktab))) tab <- cbind.data.frame(tab, ktab[[i]])
  moy <- matrix(apply(tab,2,mean),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)
  tab <- tab-moy
  ecart <- matrix(apply(tab,2,var),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)*(nrow(tab)-1)/nrow(tab)
  if (scale.unit==TRUE) tab <- tab / sqrt(ecart)
  ktab.aux <- ktab.data.frame(tab,blocks=ktab$blo,tabnames=tab.names(ktab))
  
  if (is.ktab(ktab.illu) ==TRUE){
    tab.illu <- matrix(0,nrow(ktab.illu[[1]]),0)
    for (i in 1:length(tab.names(ktab.illu))) tab.illu <- cbind.data.frame(tab.illu, ktab.illu[[i]])
    moy <- matrix(rep(moy[1,],nrow(tab.illu)),nrow=nrow(tab.illu),ncol=ncol(tab.illu),byrow=TRUE)
    tab.illu <- tab.illu-moy
    ecart.illu <- matrix(rep(ecart[1,],nrow(tab.illu)),nrow=nrow(tab.illu),ncol=ncol(tab.illu),byrow=TRUE)
    if (scale.unit==TRUE) tab.illu <- tab.illu/sqrt(ecart.illu)
  }

  colw <- NULL
  for (i in 1:length(tab.names(ktab))) colw <- c(colw,rep(1*poids[i]/eigen(cov(ktab.aux[[i]]),sym=TRUE)$values[1],ktab.aux$blo[i]))
  if (scale.unit==TRUE) colw <- colw/(nrow(tab)-1)*nrow(tab)
  if (length(tab.names(ktab))==1) colw <- colw/colw
  eig <- eigen(cov(as.matrix(tab)%*%diag(sqrt(colw))),sym=TRUE)
  scores <- as.matrix(tab)%*%diag(sqrt(colw))%*%eig$vectors
  scores.partiel <- NULL
  for (i in 1:length(tab.names(ktab))) {
    nb.prec <- nb.column(ktab,i)$prec
    nb.suiv <- nb.column(ktab,i)$suiv
    if (i==1) score.aux <- cbind(matrix(0,nrow(ktab[[1]]),nb.prec),as.matrix(tab[,1:ktab$blo[1]]),matrix(0,nrow(ktab[[i]]),nb.suiv))%*%diag(sqrt(colw))%*%eig$vectors
    if (i!=1) score.aux <- cbind(matrix(0,nrow(ktab[[i]]),nb.prec),as.matrix(tab[,(sum(ktab$blo[1:(i-1)])+1):sum(ktab$blo[1:i])]),matrix(0,nrow(ktab[[i]]),nb.suiv))%*%diag(sqrt(colw))%*%eig$vectors
#    scores.partiel <- rbind(scores.partiel,score.aux)
    scores.partiel <- rbind(scores.partiel,score.aux/poids[i]*sum(poids))
  }
#  scores.partiel <- scores.partiel * length(tab.names(ktab))
  MFASENSO <- list()
  MFASENSO$moyen <- scores[,1:nbcoord]
  MFASENSO$partiel <- matrix(0,nrow(ktab[[1]]),0)
  for (i in 1:length(tab.names(ktab))) MFASENSO$partiel <- cbind(MFASENSO$partiel,scores.partiel[((i-1)*nrow(ktab[[1]])+1):(i*nrow(ktab[[1]])),1:nbcoord])
  
  if (is.ktab(ktab.illu) ==TRUE){
    scores.illu <- as.matrix(tab.illu)%*%diag(sqrt(colw))%*%eig$vectors
    MFASENSO$moyen.illu <- NULL
    MFASENSO$partiel.illu <- NULL
    scores.partiel.illu <- NULL
    for (i in 1:length(tab.names(ktab))) {
      nb.prec <- nb.column(ktab,i)$prec
      nb.suiv <- nb.column(ktab,i)$suiv
      if (i==1) score.aux <- cbind(matrix(0,nrow(ktab.illu[[1]]),nb.prec),as.matrix(tab.illu[,1:ktab$blo[1]]),matrix(0,nrow(ktab.illu[[i]]),nb.suiv))%*%diag(sqrt(colw))%*%eig$vectors
      if (i!=1) score.aux <- cbind(matrix(0,nrow(ktab.illu[[i]]),nb.prec),as.matrix(tab.illu[,(sum(ktab$blo[1:(i-1)])+1):sum(ktab$blo[1:i])]),matrix(0,nrow(ktab.illu[[i]]),nb.suiv))%*%diag(sqrt(colw))%*%eig$vectors
#      scores.partiel.illu <- rbind(scores.partiel.illu,score.aux)
      scores.partiel.illu <- rbind(scores.partiel.illu,score.aux/poids[i]*sum(poids))
    }
#    scores.partiel.illu <- scores.partiel.illu * length(tab.names(ktab))
    MFASENSO$moyen.illu <- scores.illu[,1:nbcoord]
    MFASENSO$partiel.illu <- matrix(0,nrow(ktab.illu[[1]]),0)
    for (i in 1:length(tab.names(ktab.illu))) MFASENSO$partiel.illu <- cbind(MFASENSO$partiel.illu,scores.partiel.illu[((i-1)*nrow(ktab.illu[[1]])+1):(i*nrow(ktab.illu[[1]])),1:nbcoord])
 }
  return(MFASENSO)
}
