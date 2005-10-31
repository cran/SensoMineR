pca <- function (df,supind=NULL,supvar=NULL,row.w=rep(1, nrow(df)-length(supind))/(nrow(df)-length(supind)),scale.unit=TRUE,coord=c(1,2),graph=TRUE,main.title=NULL,clabel=1,cex=0.7,font=1,csub = 1,col="black",lty=1) {

#library(MASS)

ginv<-function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2:1])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}

coordindsup=NULL
coordvarsup=NULL

if (any(is.na(df))){
   for (j in 1:ncol(df)) df[,j] <- replace(df[,j],is.na(df[,j]),mean(df[,j],na.rm=TRUE))
                   }
if ((length(supind)==0)&(length(supvar)==0)) dfa <- df
if ((length(supind)==0)&(length(supvar)!=0)) dfa <- df[,-supvar]
if ((length(supind)!=0)&(length(supvar)==0)) dfa <- df[-supind,]
if ((length(supind)!=0)&(length(supvar)!=0)) dfa <- df[-supind,-supvar]

rank <- min(nrow(dfa)-1,ncol(dfa))

un<-rep(1,nrow(dfa))
g<-t(dfa)%*%diag(row.w)%*%un
X<-dfa-un%*%t(g)
X<-as.matrix(X)

if (scale.unit==TRUE) {
        InvV<-diag(1/sqrt(diag(t(X)%*%diag(row.w)%*%X)))
        X<-X%*%InvV   }
Y<-t(X)%*%diag(row.w)%*%X
dec<-svd(Y)

coordind<-X%*%dec$u
coordvar<-t(X)%*%diag(row.w)%*%coordind%*%sqrt(ginv(diag(dec$d)))
rownames(coordvar)=colnames(dfa)

if(length(supind) != 0) {
un<-rep(1,length(supind))
if (length(supvar)==0) Xsup<-df[supind,]-un%*%t(g)
if (length(supvar)>0) Xsup<-df[supind,-supvar]-un%*%t(g)
Xsup<-as.matrix(Xsup)
if (scale.unit==TRUE) Xsup<-Xsup%*%InvV
coordindsup<-Xsup%*%dec$u
colnames(coordindsup)=paste("Comp",1:ncol(dfa))
                        }

if(length(supvar) != 0) {
  un<-rep(1,nrow(dfa))
  if (length(supind)==0){
    gsup<-t(df[,supvar])%*%diag(row.w)%*%un
    Xsup<-df[,supvar]-un%*%t(gsup)
  }
  if (length(supind)>0){
    gsup<-t(df[-supind,supvar])%*%diag(row.w)%*%un
    Xsup<-df[-supind,supvar]-un%*%t(gsup)
  }
  Xsup<-as.matrix(Xsup)
  if (scale.unit==TRUE) {
    if (length(supvar)==1) Xsup<-Xsup/rep(sqrt(var(Xsup)*(nrow(Xsup)-1)/nrow(Xsup)),nrow(Xsup))
    if (length(supvar)>1){
      InvVsup<-diag(1/sqrt(diag(t(Xsup)%*%diag(row.w)%*%Xsup)))
      Xsup<-Xsup%*%InvVsup  
    }
  }
  coordvarsup<-t(Xsup)%*%diag(row.w)%*%coordind%*%sqrt(ginv(diag(dec$d)))
  colnames(coordvarsup)=paste("Comp",1:ncol(dfa))
}

coordvar=rbind(coordvar,coordvarsup)
if (length(supvar)>0) rownames(coordvar)=c(colnames(df)[-supvar],colnames(df)[supvar])
respcaindiv=rbind(coordind,coordindsup)
rownames(respcaindiv)=c(rownames(dfa),rownames(df)[supind])
colnames(respcaindiv)=paste("Comp",1:ncol(dfa))

if (graph){
  plot(coordind[, coord], asp=1,cex = 0.8,pch=20,xlim=c(min(respcaindiv[,coord[1]])*1.1,max(respcaindiv[,coord[1]])*1.1),ylim=c(min(respcaindiv[,coord[2]])*1.1,max(respcaindiv[,coord[2]])*1.1),xlab=paste("Comp",coord[1]," (",signif(dec$d[coord[1]]*100/sum(dec$d),4),"%)",sep=""),ylab=paste("Comp", coord[2]," (",signif(dec$d[coord[2]]*100/sum(dec$d),4),"%)",sep=""))
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  if (length(supind) != 0) points(coordindsup[, coord[1]],coordindsup[, coord[2]], cex = 0.8,pch=20,col="blue")
  text(coordind[, coord[1]], coordind[, coord[2]],labels = row.names(coordind), cex = clabel, pos = 4,offset = 0.2)
  if (length(supind) != 0) text(coordindsup[, coord[1]], coordindsup[, coord[2]],labels = row.names(coordindsup), cex = clabel, pos = 4,offset = 0.2,col="blue")
  if (length(supind) == 1) text(respcaindiv[nrow(respcaindiv), coord[1]], respcaindiv[nrow(respcaindiv), coord[2]],labels = row.names(respcaindiv)[nrow(respcaindiv)], cex = clabel, pos = 4,offset = 0.2,col="blue")
  if (length(main.title)!=0) title(main = main.title,cex.sub = cex,font.sub = font, col.sub = col, adj = 0.5)
  if (length(main.title)==0) title(main = "Individuals factor map", cex.sub = cex,font.sub = font, col.sub = col, adj = 0.5)
  get(getOption("device"))()

  if (scale.unit==TRUE){
    if (length(supvar)==0) senso.corcircle(coordvar[, coord])
    if (length(supvar)>0){
      senso.corcircle(coordvar[1:(nrow(coordvar)-length(supvar)), coord])
      if (length(supvar) == 1) senso.corcircle(matrix(rbind(coordvar[nrow(coordvar), coord],coordvar[nrow(coordvar), coord]),nrow=2,dimnames=list(c(" ",row.names(coordvar)[nrow(coordvar)]),NULL)),col="blue",add=TRUE,lty=2)
      if (length(supvar) > 1) senso.corcircle(coordvar[(nrow(coordvar)-length(supvar)+1):nrow(coordvar), coord],col="blue",lty=2,add=TRUE)
    }
    aaa = paste("Correlation circle: ", "comp",coord[1]," (",signif(dec$d[coord[1]]*100/sum(dec$d),4),"%) - ", "comp", coord[2]," (",signif(dec$d[coord[2]]*100/sum(dec$d),4),"%)",sep="")
    if (length(main.title)!=0) title(main = main.title,sub=aaa, , cex.sub = cex,font.sub = font, col.sub = col,adj=0.5,line=2.8)
    if (length(main.title)==0) title(main = aaa, cex.sub = cex,font.sub = font, col.sub = col,adj=0.5,line=2.8)
  }
  else {
    if (length(supvar)==0) senso.arrow(coordvar[, coord],lab = rownames(coordvar))
    if (length(supvar)>0){
      senso.arrow(coordvar[1:(nrow(coordvar)-length(supvar)), coord],lab = rownames(coordvar)[1:(nrow(coordvar)-length(supvar))])
      if (length(supvar) == 1) senso.arrow(matrix(rbind(coordvar[nrow(coordvar), coord],coordvar[nrow(coordvar), coord]),nrow=2,dimnames=list(c(" ",row.names(coordvar)[nrow(coordvar)]),NULL)),col="blue",lty=2,add.plot=TRUE)
      if (length(supvar) > 1) senso.arrow(coordvar[(nrow(coordvar)-length(supvar)+1):nrow(coordvar), coord],col="blue",lty=2,add.plot=TRUE)
    }
    aaa=paste("Variables factor map: ", "comp",coord[1]," (",signif(dec$d[coord[1]]*100/sum(dec$d),4),"%) - ", "comp", coord[2]," (",signif(dec$d[coord[2]]*100/sum(dec$d),4),"%)",sep="")
    if (length(main.title)!=0) title(main = main.title,sub=aaa, cex.sub = cex,font.sub = font, col.sub = col,adj=0.5,line=1.8)
    if (length(main.title)==0) title(main = aaa, cex.sub = cex,font.sub = font, col.sub = col,adj=0.5,line=1.8)
  }
}

colnames(coordind)=colnames(coordvar)=paste("Comp",1:ncol(dfa))
resacp = list()
resacp$eig = dec$d[1:rank]
resacp$li = coordind[,1:rank]
resacp$co = coordvar[1:(nrow(coordvar)-length(supvar)),1:rank]
if (length(supind)==1) resacp$lisup = matrix(coordindsup[,1:rank],nrow=1,dimnames=list(rownames(coordindsup),colnames(coordindsup)[1:rank]))
if (length(supind)>1) resacp$lisup = coordindsup[,1:rank]
if (length(supvar)==1) resacp$cosup = matrix(coordvar[nrow(coordvar),1:rank],nrow=1,dimnames=list(rownames(coordvar)[nrow(coordvar)],colnames(coordvar)[1:rank]))
if (length(supvar)>1) resacp$cosup = coordvar[(nrow(coordvar)-length(supvar)+1):nrow(coordvar),1:rank]
return(resacp)

}
