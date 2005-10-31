"paneliperf" <- function(donnee,formul,formul.j="~Product",col.j,firstvar,lastvar=ncol(donnee),synthesis=FALSE,random=TRUE,graph=TRUE){

  for (j in 1 :(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])
  formul.j = as.formula(formul.j)
  lab <- labels(donnee)[[2]]
  equation <- as.character(formul.j)
  lab.j <- levels(donnee[,col.j])

  dim.donnee <- ncol(donnee)
  for (i in 1:dim.donnee) {
    if (gsub(" ","",strsplit(equation,split="+",extended=FALSE)[[2]][1])==lab[i]) col.p <- i
  }
  r2 <- prob <- vtest <- res <- matrix(0,length(levels(donnee[,col.j])),lastvar+1-firstvar)

  for (j in 1:length(lab.j)){
    for (varendo in firstvar:lastvar) {
      formule <- paste(lab[varendo],"~ C(")
      aux2 <- equation[2]
      aux3 <- strsplit(aux2,"+",extended=FALSE)[[1]]

      for (i in 1:length(aux3)) {
        if (any(grep("%",aux3[i]))) {
              formule <- paste(formule,strsplit(aux3[i],"%",extended=FALSE)[[1]][1],",sum)%in%C(",strsplit(aux3[i],"%",extended=FALSE)[[1]][3],",sum)")
        }
        else
        if (any(grep(":",aux3[i]))) {
           formule <- paste(formule,strsplit(aux3[i],":",extended=FALSE)[[1]][1],",sum) : C(",strsplit(aux3[i],":",extended=FALSE)[[1]][2],",sum)")
        }
        else {
        formule <- paste(formule,aux3[i],",sum)")
        }
        if (i < length(aux3))  formule <- paste (formule, "+C(")
      }

      formule <- as.formula(formule)
      aux <- summary(aov( formule , data = donnee, subset=(donnee[,col.j]==levels(donnee[,col.j])[j]),na.action =na.omit))[[1]]
      prob[j,varendo-firstvar+1] <- aux[1,5]
      vtest[j,varendo-firstvar+1] <- -qnorm(aux[1,5]/2)
      res[j,varendo-firstvar+1] <- sqrt(aux[nrow(aux),3])
      r2[j,varendo-firstvar+1] <- summary.lm(aov( formule , data = donnee, subset=(donnee[,col.j]==levels(donnee[,col.j])[j]),na.action =na.omit))$r.squared

    }
}

  dimnames(prob) <-   dimnames(r2) <-    dimnames(vtest) <-  dimnames(res) <- list(levels(donnee[,col.j]),labels(donnee)[[2]][firstvar:lastvar])

  aa <- averagetable(donnee,formul=formul,firstvar=firstvar)
  agree <- matrix(0,length(lab.j),lastvar-firstvar+1)
  for (j in 1:length(lab.j)){
      bb <- averagetable(donnee,formul=formul.j,subset=(donnee[,col.j]==lab.j[j]),firstvar=firstvar,lastvar=lastvar)
      agree[j,] <- diag(cor(aa,bb))
    }
  dimnames(agree) <- list(lab.j,lab[firstvar:lastvar])  

 if (graph){
  aux.agree=agree
  colnames(aux.agree)=paste(colnames(agree),".a",sep="")
  aux.vtest=-qnorm(prob/2)
  colnames(aux.vtest)=paste(colnames(prob),".p",sep="")
  afmult(ktab.data.frame(cbind.data.frame(aux.agree,aux.vtest),blocks=c(ncol(aux.agree),ncol(aux.vtest)),tabnames=c("Agree","Prob")),scannf=FALSE,option="lambda1")
 }
  paneliperf = list() 
  paneliperf$prob.ind = prob
  paneliperf$vtest.ind = vtest
  paneliperf$res.ind = res
  paneliperf$r2.ind = r2
  paneliperf$signif.ind=apply(matrix(as.integer(prob<0.05),nrow=nrow(prob)),1,sum,na.rm=TRUE) # nb de descripteur significatif / juge
  paneliperf$agree.ind = agree
  
  if (synthesis){
    bb=panelperf(donnee,formul=formul,subset=NULL,firstvar=firstvar,lastvar=lastvar,random=random)
    aux=cbind(-qnorm(bb$p.value/2),apply(agree,2,median,na.rm=T),apply(prob,2,median,na.rm=T),bb$res,bb$r2)
    colnames(aux)=c(colnames(bb$p.value),"median(agree)","median(prob.ind)","stdev Res","R2")
    if (graph){
      get(getOption("device"))()
      pca(aux,supvar=(ncol(aux)-3):ncol(aux),main.title="PCA on the P-values issued from the AOV model")
    }
    paneliperf$complete = aux
    paneliperf$p.value = bb$p.value
    paneliperf$variability = bb$variab
    paneliperf$res = bb$res
  }
  return(paneliperf)
}
