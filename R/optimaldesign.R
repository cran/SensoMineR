optimaldesign <- function(nbPanelist,nbProd,nbProdByPanelist=nbProd, seed=NULL){

  if (is.null(seed)) seed <- sample(1:1000,1)
  set.seed(seed)

  mat <- data.frame(prod = as.factor(rep(1:nbProd,nbPanelist)),juge = as.factor(rep(1:nbPanelist,each=nbProd)))
  desD <- optFederov(~prod+juge,mat,nTrials=nbPanelist*nbProdByPanelist,maxIteration=1000,nRepeats=1000)$design

  comp <- data.frame(prod=as.factor(rep(desD[,1],nbProdByPanelist)),juge=as.factor(rep(desD[,2],nbProdByPanelist)),rang=as.factor(rep(1:nbProdByPanelist,each=nrow(desD))))
  comp=comp[sample(1:nrow(comp)),]
  desF <- optFederov(~juge+rang+prod,comp,nTrials=nbPanelist*nbProdByPanelist,maxIteration=1000,nRepeats=1000)$design
  desF <- desF[order(desF[,2],desF[,3]),]

  design <- matrix(NA,nbPanelist,nbProdByPanelist)
  for (j in 1:nbPanelist){
    if (setequal(desF[desF[,2]==j,3],1:nbProdByPanelist)) design[j,] <- desF[desF[,2]==j,1]
    else {
      design[j,desF[desF[,2]==j,3]] <- desF[desF[,2]==j,1]
	  design[j,!((1:nbProdByPanelist)%in%desF[desF[,2]==j,3])] <- desF[desF[,2]==j,][which(duplicated(desF[desF[,2]==j,3],fromLast=TRUE)),1]
    }
  }
  rownames(design) = paste("Panelist",1:nbPanelist)
  colnames(design) = paste("Rank",1:nbProdByPanelist)

  valid.succ = matrix(0,nbProd,nbProd)
  for (i in 1:nrow(design)){
    for (k in 1:(nbProdByPanelist-1)){
     valid.succ[design[i,k],design[i,k+1]]=valid.succ[design[i,k],design[i,k+1]]+1
    }
  }

  tab = data.frame(prod=as.factor(design),juge=as.factor(rep(1:nbPanelist,nbProdByPanelist)),rank=as.factor(rep(1:nbProdByPanelist,each=nbPanelist)))
  res = list(design=design, valid.rank=table(tab[,1],tab[,3],dnn=list("Product","Rank")), valid.succ=valid.succ, ProdPanelist = table(tab[,1],tab[,2],dnn=list("Product","Panelist")), RankPanelist = table(tab[,2],tab[,3],dnn=list("Panelist","Rank")))
  return(res)
}