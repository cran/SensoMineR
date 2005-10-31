"construct.axes" <- function(ktableau,coord=c(1,2),scale.unit=TRUE,centerbypanelist=FALSE,scalebypanelist=FALSE){

  nbcoord=max(coord)
  matrice <- ktableau[[1]]
  matrice[,1]=as.factor(matrice[,1])
  matrice[,2]=as.factor(matrice[,2])
  for (i in 2:length(ktableau$blo))  matrice <- cbind.data.frame(matrice,ktableau[[i]])

  oo <- order(matrice[,2])
  matrice <- matrice[oo,]
  oo <- order(matrice[,1])
  matrice <- matrice[oo,]

  nbjuge <- sum(as.integer(summary(matrice[,1])!=0))
  nbprod <- length(levels(matrice[,2]))
  nbdesc <- dim(matrice)[2]-2

  moy.aux=scalebypanelist(matrice,col.j=1,col.p=2,firstvar=3,center=centerbypanelist,scale=scalebypanelist)
  
  ###AF with active data the averages for all the panelist 
  nbactif <- nbprod
  ktable <- ktab.data.frame(moy.aux[1:nbactif,-c(1,2)],blocks=ktableau$blo[-1],tabnames=tab.names(ktableau)[-1])
  ktable2 <- ktab.data.frame(moy.aux[(nbactif+1):dim(moy.aux)[1],-c(1,2)],blocks=ktableau$blo[-1],tabnames=tab.names(ktableau)[-1])
  res.afm <- mfasenso(ktable,ktable2,scale.unit=scale.unit,nbcoord=nbcoord)

  axe <- list()
  axe$moyen <- data.frame(rbind(res.afm$moyen,res.afm$moyen.illu),as.factor(moy.aux[,2]),as.factor(moy.aux[,1]))
  dimnames(axe$moyen)[2][[1]]<-c (paste("Comp", 1:nbcoord, sep = ""),"Product","Panelist")

  if (length(ktableau$blo)>2) {
    axe$partiel <- data.frame(rbind(res.afm$partiel,res.afm$partiel.illu),as.factor(moy.aux[,2]),as.factor(moy.aux[,1]))
    dimnames(axe$partiel)[2][[1]][(dim(axe$partiel)[2]-1):dim(axe$partiel)[2]] <- c("Product","Panelist")
    for (i in 1:length(ktable$blo)) dimnames(axe$partiel)[2][[1]][((i-1)*nbcoord+1):(i*nbcoord)]<-paste("Comp", 1:nbcoord, sep = "",tab.names(ktable)[i])
  }
  aa=cor(moy.aux[1:nbprod,-(1:2)],axe$moyen[1:nbprod,coord])
  get(getOption("device"))()
  senso.corcircle(aa, fullcircle=TRUE)
  title(main = paste("Correlation circle (comp ",coord[1]," - comp ",coord[2],")",sep=""))
  get(getOption("device"))()
  par(mar = c(0,0,2,0))
  senso.label( axe$moyen[1:nbprod,coord], clabel=0, cpoint=0.8, include.origin = FALSE)
  text( axe$moyen[1:nbprod,coord[1]], axe$moyen[1:nbprod,coord[2]], labels = levels( axe$moyen[1:nbprod,nbcoord+1]), cex = 0.8, pos = 4, offset = 0.2)
  title(main = paste("Individuals factor map (comp ",coord[1]," - comp",coord[2],")",sep=""))

  return(axe)
}
