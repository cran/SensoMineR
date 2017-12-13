JAR <- function(x, col.p, col.j, col.pref, jarlevel="jar"){

  fct.delete.first <- function(x){x[-1]}
  ind.jar <- (1:ncol(x))[-c(col.p,col.j,col.pref)]
  nbjar <- length(ind.jar)
  for (i in ind.jar){
    x[,i]=relevel(x[,i],jarlevel)
  }

  nbmod=rep(0,ncol(x))
  for(j in ind.jar){
    nbmod[j]=nlevels(x[,j])-1
  }
  nbmodtot=sum(nbmod)
  varesid=rep(1,nbjar+1)

  penal1 <- penal2 <- matrix(NA,nrow=nbmodtot,ncol=3)
  nommod=rep("a",nbmodtot)
  colnames(penal1)=c("One-dimension penalty (all products)","Standard error","p-value")

  ifin=0
  for(j in ind.jar){
    res=lm(x[,col.pref]~x[,j])
    varesid[j]= anova(res)[2,3]
    ideb=ifin+1
    ifin=ideb+nbmod[j]-1
    npar=nbmod[j]+1
    penal1[c(ideb:ifin),1]=-res$coefficients[2:npar]
    penal1[c(ideb:ifin),2]=summary(res)[[4]][2:npar,2]
    nommod[c(ideb:ifin)]=levels(x[,j])[2:npar]
    penal1[c(ideb:ifin),3]=summary(res)[[4]][2:npar,4]
  }
 res <- lm(x[,col.pref]~.,data=x[,c(ind.jar,col.p,col.j)])
 penal2 <- summary(res)$coef[2:(nbmodtot+1),c(1,2,4)]
 colnames(penal2)[1]="Multidimensional penalty (all products)"
 rownames(penal1)=rownames(penal2)=nommod
 Frequency <- matrix(NA,nrow=nbmodtot,ncol=nlevels(x[,col.p]))
 for (j in 1:ncol(Frequency)) Frequency[,j] <- unlist(lapply(lapply(x[x[,col.p]==levels(x[,col.p])[j],ind.jar],table),fct.delete.first))
 Frequency <- sweep(Frequency,2,table(x[,col.p]),FUN="/")*100
 rownames(Frequency)=nommod
 colnames(Frequency)=levels(x[,col.p])

 res <- list(penalty1 = penal1,penalty2 = penal2, Frequency=Frequency)
 class(res) <- c("JAR", "list ")
 return(res)
}
