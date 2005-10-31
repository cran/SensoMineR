"panellipse.session" <- function(donnee,col.p,col.j,col.s,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,level.search.desc=0.5,centerbypanelist=TRUE,scalebypanelist=FALSE){

for (j in 1:(firstvar-1)) donnee[,j]=as.factor(donnee[,j])
labseance=levels(donnee[,col.s])
nbseance <- length(labseance)
nbprod <- length(levels(donnee[,col.p]))
nbjuge <- length(levels(donnee[,col.j]))

donnee <- search.desc(donnee,col.j=col.j,col.p=col.p,firstvar=firstvar,lastvar=lastvar,level=level.search.desc)
if (nbseance <2) print("This procedure is not adapted, there is only one session")
  oo=order(donnee[,col.j])
  donnee<- donnee[oo,]
  oo=order(donnee[,col.p])
  donnee<- donnee[oo,]
  oo=order(donnee[,col.s])
  donnee<- donnee[oo,]

don <- cbind.data.frame(donnee[donnee[,col.s]==labseance[1],col.j],donnee[donnee[,col.s]==labseance[1],col.p])

for (seance in 1:nbseance)  don <- cbind.data.frame(don,data.frame(donnee[donnee[,col.s]==labseance[seance],firstvar:ncol(donnee)],row.names=paste(donnee[donnee[,col.s]==labseance[seance],col.p],donnee[donnee[,col.s]==labseance[seance],col.j],sep=".")))
  colnames(don) <- colnames(donnee)[c(col.j,col.p,rep(firstvar:ncol(donnee),nbseance))]
  colnames(don) <- paste(colnames(don),c("","",rep(paste(".S",1:nbseance,sep=""),rep(ncol(donnee)-firstvar+1,nbseance))),sep="")
  aa <- matrix(0,ncol(donnee)-firstvar+1,2)
  res.average=averagetable(don,formul=as.formula(paste("~",colnames(don)[2])),firstvar=3)
  for (j in 1:(ncol(donnee)-firstvar+1)) aa[j,1]<-pca(res.average[,(ncol(donnee)-firstvar+1)*(0:(nbseance-1))+j],graph=FALSE)$eig[1]
  aa[,2]=aa[,1]/nbseance*100
  rownames(aa) <- colnames(donnee[,firstvar:ncol(donnee)])
  colnames(aa) <- c("eig1","Reproductibility")
  aa=aa[rev(order(aa[,1])),]
  ktab.donnee <- ktab.data.frame(don,blocks=c(2,rep(ncol(donnee)-firstvar+1,nbseance)),tabnames=c("JP",paste("S",1:nbseance,sep="")))
  axe <- construct.axes(ktab.donnee,coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist)
  plotpanelist(axe$moyen,coord=coord)
  get(getOption("device"))()
  simul <- simulation(axe,nbbloc=length(ktab.donnee$blo)-1,nbchoix=nbchoix,nbsimul=nbsimul)
  plotellipse(simul,alpha=alpha,coord=coord)  
  return(aa)
}
