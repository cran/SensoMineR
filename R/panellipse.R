"panellipse" <- function(donnee,col.p,col.j,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,bloc=NULL,name.bloc=NULL,level.search.desc=0.5,centerbypanelist=TRUE,scalebypanelist=FALSE){
  don <- cbind.data.frame(donnee[,col.j],donnee[,col.p],donnee[,firstvar:lastvar])
  colnames(don) <- colnames(donnee)[c(col.j,col.p,firstvar:lastvar)]
  if (length(bloc)<2) ktab.donnee <- ktab.data.frame(don,blocks=c(2,lastvar-firstvar+1),tabnames=c("JP","Gr1"))
  if (length(bloc)>1) ktab.donnee <- ktab.data.frame(don,blocks=c(2,bloc),tabnames=c("JP",name.bloc))
  ktab.interesting.desc <- search.desc.ktab(ktab.donnee,level=level.search.desc)
  axe <- construct.axes(ktab.interesting.desc,coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist)
  plotpanelist(axe$moyen,coord=coord)
  get(getOption("device"))()
  simul <- simulation(axe,nbbloc=length(ktab.interesting.desc$blo)-1,nbchoix=nbchoix,nbsimul=nbsimul)
  plotellipse(simul,alpha=alpha,coord=coord)  
}
