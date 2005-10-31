"panelmatch" <- function(donnee,col.p,col.j,firstvar,alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,bloc=NULL,name.bloc=NULL,centerbypanelist=TRUE,scalebypanelist=FALSE,name.panelist=FALSE,cex=1,color=NULL,hierar=NULL){

hotelling <- function(d1,d2,n1=nrow(d1),n2=nrow(d2)){
    k <- ncol(d1)
    xbar1 <- apply(d1,2,mean)
    xbar2 <- apply(d2,2,mean)
    dbar <- xbar2-xbar1
    v <- ((n1-1)*var(d1)+(n2-1)*var(d2))/(n1+n2-2)
    t2 <- n1*n2*dbar%*%solve(v)%*%dbar/(n1+n2)
    f <- (n1+n2-k-1)*t2/((n1+n2-2)*k)
    return(pf(f,k,n1+n2-k-1,lower.tail=FALSE))
}


htabdes <- function(X) {
    nbnivh <- length(X)
    nbvarh <- X
    if (nbnivh>1) {
    for (i in 2:nbnivh) {
        for (j in 1:length(X[[i]])) {    
        nbvarh[[i]][j] <- 0
        if (j==1)   for (k in 1:X[[i]][1])    nbvarh[[i]][j] <- nbvarh[[i]][j]+nbvarh[[i-1]][k]
        else {         
            a <- 0     
        b <- 0      
            for (n in 1:(j-1))   a <- a+X[[i]][n]
        a <- a+1
        for (n in 1:j)   b <- b+X[[i]][n]
            for (k in a:b)   nbvarh[[i]][j] <- nbvarh[[i]][j]+nbvarh[[i-1]][k]
            }       
        }
    }
                  }
    if (nbnivh==1) nbvarh=X
    return(nbvarh)
}

hdil <- function(X) {
    nbnivh <- length(X)
    dil <- X
    a <- NULL
    dil[[nbnivh]] <- rep(length(X[[nbnivh]]),length(X[[nbnivh]]))
    if (nbnivh>1) {
    for (i in 1:(nbnivh-1)) { 
    h <- nbnivh-i
    k <- nbnivh-i+1
        for (j in 1:length(X[[k]])) {
            a <- c(a,rep( X[[k]][j]*dil[[k]][j],X[[k]][j] ))
        }
        dil[[h]] <- a
        a <- NULL
    }
                  }
    return(dil)
}


hweight <- function (X, H) {
    nbnivh <- length(H)
    nivo <- htabdes(H)
    cw <- rep(1, ncol(X))
    cw.partiel <- H
    for (n in 1:nbnivh) {
    Xinter <- ktab.data.frame(X,nivo[[n]],w.col=cw)
        sepan <- sepan(Xinter)
        nbloc <- length(sepan$blo)
    rank.fac <- factor(rep(1:nbloc, sepan$rank))
        cwinter <- NULL
    for (i in 1:nbloc) {
            cwinter <- c(cwinter, rep(1/sepan$Eig[rank.fac == i][1],sepan$blo[i]))
    }
        cw <- cw * cwinter
        cw.partiel[[n]] <- cw
    }
    return(cw.partiel)
}

  if (length(color)==0) color = c("black","red","green3","blue",
    "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
    "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
    "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
    "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")

 nvar=NULL
 nvarpanel=NULL
 namevar=NULL
 namegr=NULL
 aux2=NULL
 aamoy=NULL
 numligne=NULL
 nbjuge=NULL
 nbpanel=length(donnee)
 labprod = levels(as.factor(donnee[[1]][,col.p]))
 nbprod=length(labprod)
 for (bloc in 1:length(donnee)) {
   nvar=c(nvar,ncol(donnee[[bloc]])-firstvar+1)
 }
 if (length(hierar)!=0) nvar=hierar[[1]]
 if (length(hierar)==0) hierar=list(nvar)
 for (bloc in 1:nbpanel){
#    aux <- search.desc(donnee[[bloc]],level=level.search.desc,col.j=col.j,col.p=col.p,firstvar=firstvar) # pas de sélection car sinon la hiérarchie est à revoir
    aux = scalebypanelist(donnee[[bloc]],center=centerbypanelist,scale=scalebypanelist,col.j=col.j,col.p=col.p,firstvar=firstvar)
    nbjuge=c(nbjuge,nrow(aux)/nbprod-1)
    numligne=c(numligne,nrow(aux)-nbprod)
    namegr=c(namegr,paste("P",bloc,sep=""))
    namevar=c(namevar,paste("P",bloc,"_",colnames(donnee[[bloc]])[-(1:(firstvar-1))],sep=""))
    nvarpanel=c(nvarpanel,ncol(donnee[[bloc]])-firstvar+1)
    if (bloc==1) aamoy = aux[1:nbprod,]
    if (bloc!=1) aamoy = cbind(aamoy,aux[1:nbprod,-(1:2)])
 }
 rownames(aamoy)=aamoy[,2]
 tab=aamoy[,-(1:2)]
 colnames(tab)=namevar
 moy <- matrix(apply(tab,2,mean),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)
 tab <- tab-moy
 ecart <- matrix(apply(tab,2,var),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)*(nrow(tab)-1)/nrow(tab)
 if (scale.unit==TRUE) tab <- tab / sqrt(ecart)
 aamoy = tab
 resafmh=hmfa(as.data.frame(aamoy),hierar,coord=coord)
 dilat = hdil(hierar)

##
for (j in 1:length(hierar)){
  if (length(htabdes(hierar)[[j]]) == nbpanel) nivpanel = j
}
poids = hweight(aamoy,hierar)

for (bloc in 1:nbpanel){
#    aux <- search.desc(donnee[[bloc]],level=level.search.desc,col.j=col.j,col.p=col.p,firstvar=firstvar) # pas de sélection car sinon la hiérarchie est à revoir
    aux = scalebypanelist(donnee[[bloc]],center=centerbypanelist,scale=scalebypanelist,col.j=col.j,col.p=col.p,firstvar=firstvar)
    tab=aux[1:nbprod,-(1:2)]
    moy <- matrix(apply(tab,2,mean),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)
    ecart <- matrix(apply(tab,2,var),nrow=nrow(tab),ncol=ncol(tab),byrow=TRUE)*(nrow(tab)-1)/nrow(tab)
    aux=aux[-(1:nbprod),]
    aux2 = cbind.data.frame(paste(namegr[bloc],aux[,1],sep="_"),aux[,2])
    aux = aux[,-(1:2)]
    tab.illu = aux
    moy <- matrix(rep(moy[1,],nrow(tab.illu)),nrow=nrow(tab.illu),ncol=ncol(tab.illu),byrow=TRUE)
    tab.illu <- tab.illu-moy
    ecart.illu <- matrix(rep(ecart[1,],nrow(tab.illu)),nrow=nrow(tab.illu),ncol=ncol(tab.illu),byrow=TRUE)
    if (scale.unit==TRUE) tab.illu <- tab.illu/sqrt(ecart.illu)
    aux = tab.illu
    names(aux)=NULL
    if (bloc == 1) aux3=cbind.data.frame(aux,matrix(0,nrow(aux),sum(nvarpanel[(bloc+1):length(nvarpanel)])))
    if ((bloc >1) & (bloc != nbpanel)) aux3 = data.frame(matrix(0,nrow(aux),sum(nvarpanel[1:(bloc-1)])),aux,matrix(0,nrow(aux),sum(nvarpanel[(bloc+1):length(nvarpanel)])))
    if (bloc==nbpanel) aux3 = cbind.data.frame(matrix(0,nrow(aux),sum(nvarpanel[1:(bloc-1)])),aux)
    colnames(aux3)=namevar
#    formule <- as.matrix(aux3) %*% diag(poids[[nivpanel]]) %*% t(aamoy) %*% diag(rep((1/nbprod),nbprod))
    formule <- as.matrix(aux3) %*% diag(poids[[length(hierar)]]) %*% t(aamoy) %*% diag(rep((1/nbprod),nbprod))
    formule <- formule %*% as.matrix(resafmh$hmfa$li)[,1:max(coord)] %*% diag(1/resafmh$hmfa$eig[1:max(coord)]) * dilat[[nivpanel]][bloc]
    indivbloc = cbind.data.frame(formule,aux2[,c(2,1)])
    if (bloc !=1) indiv = rbind(indiv,indivbloc)
    if (bloc ==1) indiv = indivbloc
}
    colnames(indiv)[(ncol(indiv)-1):ncol(indiv)]=c("Product","Panelist")

 comp.simul=list()
 for (bloc in 1:nbpanel){
   auxi = cbind.data.frame(indiv[1:nbprod,-((ncol(indiv)-1):ncol(indiv))],rownames(aamoy),rep(0,nbprod))
   colnames(auxi)[(ncol(indiv)-1):ncol(indiv)]=c("Product","Panelist")
   aa1=list()
   aa1$moyen=rbind.data.frame(auxi,indiv[(sum(numligne[1:bloc])-numligne[bloc] +1):(sum(numligne[1:bloc])),])
   aa1$moyen[,ncol(aa1$moyen)]=as.factor(aa1$moyen[,ncol(aa1$moyen)])
   simul <- simulation(aa1, nbbloc = 1,nbsimul=nbsimul,nbchoix=nbchoix)
   if (bloc==1) auxil = simul$moy$simul[,-ncol(simul$moy$simul)]
   if (bloc!=1) auxil=cbind(auxil,simul$moy$simul[,-ncol(simul$moy$simul)])
   if (bloc==1) comp.simul$partiel$simul=cbind.data.frame(simul$moy$simul[,-ncol(simul$moy$simul)],as.factor(paste("P",bloc,"_",simul$moy$simul[,ncol(simul$moy$simul)],sep="")))
   if (bloc!=1) comp.simul$partiel$simul=rbind.data.frame(comp.simul$partiel$simul, cbind.data.frame(simul$moy$simul[,-ncol(simul$moy$simul)],as.factor(paste("P",bloc,"_",simul$moy$simul[,ncol(simul$moy$simul)],sep=""))))
 }
 
  eig =   cbind(resafmh$hmfa$eig,signif(resafmh$hmfa$eig/sum(resafmh$hmfa$eig)*100,4))
  auxi= cbind.data.frame(resafmh$hmfa$li[,1:max(coord)],as.factor(indiv[1:nbprod,ncol(indiv)-1]),as.factor(rep(0,nbprod)))
  colnames(auxi) = colnames(indiv)
  plotpanelist(rbind.data.frame(auxi,indiv),eig=eig,color=color,coord=coord,name=name.panelist,cex=cex)
 comp.simul$moy$simul=auxil[,1:max(coord)]
 aux4=resafmh$coord.partial[[nivpanel]][[1]]
 for (bloc in 2:nbpanel) aux4=rbind(aux4,resafmh$coord.partial[[nivpanel]][[bloc]])
 for (ncoord in 1:max(coord))  comp.simul$moy$simul[,ncoord] = apply(auxil[,(max(coord)*(0:(nbpanel-1)))+ncoord],1,mean)
 comp.simul$moy$simul = cbind.data.frame(comp.simul$moy$simul,as.factor(simul$moy$simul[,ncol(comp.simul$partiel$simul)]))
 comp.simul$moy$P = cbind.data.frame(resafmh$hmfa$li[,1:max(coord)],rownames(resafmh$hmfa$li))
 comp.simul$partiel$P = cbind.data.frame(aux4[,1:max(coord)],rownames(aux4))

  if (nivpanel == length(hierar)){
  get(getOption("device"))(12,8)
  plotellipse(comp.simul,alpha=alpha,eig=eig,coord=coord,color=color,cex=cex)  
  if (length(names(donnee))==0) legend("bottomleft",legend=paste("Panel",1:nbpanel,sep=" "),lty=1:nbpanel,cex=0.8,bg="white")
  if (length(names(donnee))!=0) legend("bottomleft",legend=names(donnee),lty=1:nbpanel,cex=0.8,bg="white")
}

  if (nivpanel < length(hierar)){

    aux4=resafmh$coord.partial[[nivpanel+1]][[1]]
    for (k in 2:length(resafmh$coord.partial[[nivpanel+1]])) aux4=rbind(aux4,resafmh$coord.partial[[nivpanel+1]][[k]])
    comp.simul2 = list()
    comp.simul2$moy$P=comp.simul$moy$P
    comp.simul2$moy$simul=comp.simul$moy$simul
    comp.simul2$partiel$P = cbind.data.frame(aux4[,1:max(coord)],rownames(aux4))
    for (k in 1:length(hierar[[nivpanel+1]])) {
      aux5 = array(NA,dim=c(nbprod*nbsimul,max(coord),hierar[[nivpanel+1]][k]))
      if (k==1) lim1 = 1
      if (k!=1) lim1 = sum(hierar[[nivpanel+1]][1:(k-1)])+1
      lim2 = sum(hierar[[nivpanel+1]][1:k])
      for (j in lim1:lim2) aux5[,,j-lim1+1] = as.matrix(comp.simul$partiel$simul[(nbsimul*nbprod*(j-1)+1):(nbsimul*nbprod*j),1:max(coord)])
      aux6 = apply(aux5,1:2,mean)
      if (k==1) comp.simul2$partiel$simul=cbind.data.frame(aux6,as.factor(paste("R",k,"_",rep(labprod,each=nbsimul),sep="")))
      if (k!=1) comp.simul2$partiel$simul=rbind.data.frame(comp.simul2$partiel$simul,cbind.data.frame(aux6,as.factor(paste("R",k,"_",rep(labprod,each=nbsimul),sep=""))))
    }

    get(getOption("device"))(12,8)
    plotellipse(comp.simul2,alpha=alpha,eig=eig,coord=coord,color=color,cex=cex) 
    title(sub="Grouping panels") 
    legend("bottomleft",legend=paste("Grouping panels",1:length(hierar[[nivpanel+1]]),sep=" "),lty=1:length(hierar[[nivpanel+1]]),cex=0.8,bg="white")

    get(getOption("device"))(12,8)
    plotellipseinterhmfa(comp.simul,alpha=alpha,coord=coord,nbbloc=nbpanel,eig=eig,color=color,cex=cex,hmfa=list(hierar,resafmh$coord.partial))
    title(sub="Panels") 
    if (length(names(donnee))==0) legend("bottomleft",legend=paste("Panel",1:nbpanel,sep=" "),lty=1:nbpanel,cex=0.8,bg="white")
    if (length(names(donnee))!=0) legend("bottomleft",legend=names(donnee),lty=1:nbpanel,cex=0.8,bg="white")
  }
  mat = array(NA,dim=c(nbprod,nbprod,nbpanel+1))
  for (bloc in 1:(nbpanel+1)){
    if (bloc == (nbpanel+1)) aa = indiv
    if (bloc !=(nbpanel+1)) aa = indiv[(nbprod*(sum(nbjuge[1:bloc])-nbjuge[bloc])):(nbprod*sum(nbjuge[1:bloc])),]
    for (i in 1:(nbprod-1)){
      for (j in i:nbprod){
        if (length(nbchoix)==0){
           if (bloc <(nbpanel+1)) mat[i,j,bloc] = mat[j,i,bloc] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbjuge[bloc],nbjuge[bloc])
           if (bloc ==(nbpanel+1)) mat[i,j,bloc] = mat[j,i,bloc] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],sum(nbjuge),sum(nbjuge))
        }
        if (length(nbchoix)!=0) mat[i,j,bloc] = mat[j,i,bloc] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
    }}
    for (i in 1:nbprod) mat[i,i,bloc]=1
  }
  if (length(names(donnee))==0) dimnames(mat)=list(labprod,labprod,c(paste("P",1:nbpanel,sep=""),"global"))
  if (length(names(donnee))!=0) dimnames(mat)=list(labprod,labprod,c(names(donnee),"global"))

  mat2 = array(NA,dim=c(nbpanel,nbpanel,nbprod))
  for (i in 1:nbprod){
    for (k in 1:(nbpanel-1)){
      aa = indiv[(nbprod*(sum(nbjuge[1:k])-nbjuge[k])):(nbprod*sum(nbjuge[1:k])),]
      for (kk in (k+1):nbpanel){
        aa2 = indiv[(nbprod*(sum(nbjuge[1:kk])-nbjuge[kk])):(nbprod*sum(nbjuge[1:kk])),]
        if (length(nbchoix)==0) mat2[k,kk,i] = mat2[kk,k,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa2[aa2[,ncol(aa2)-1]==labprod[i],coord])
        if (length(nbchoix)!=0) mat2[k,kk,i] = mat2[kk,k,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa2[aa2[,ncol(aa2)-1]==labprod[i],coord],nbchoix,nbchoix)
      }}
    for (k in 1:nbpanel) mat2[k,k,i]=1
  }
  if (length(names(donnee))==0) dimnames(mat2)=list(paste("Panel",1:nbpanel,sep=" "),paste("Panel",1:nbpanel,sep=" "),labprod)
  if (length(names(donnee))!=0) dimnames(mat2)=list(names(donnee),names(donnee),labprod)

  res <- resafmh
  res$hotelling=list(bypanel=mat,byproduct=mat2)
  res$indiv=indiv
  return(res)
}
