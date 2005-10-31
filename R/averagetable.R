"averagetable" <- function(donnee,formul,subset=NULL,method="coeff",firstvar,lastvar=ncol(donnee),file=NULL){
 
if ((method!="coeff")&(method!="mean")) stop(paste("The method",method,"is unknown. Use coeff or mean"))

for (j in 1:(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])

formul = as.formula(formul)
lab <- labels(donnee)[[2]]
equation <- as.character(formul)
dim.donnee <- ncol(donnee)
for (i in 1:dim.donnee) {
  if (gsub(" ","",strsplit(equation,split="+",extended=FALSE)[[2]][1])==lab[i]) col.p <- i
}
nbprod <- length(levels(donnee[,col.p]))
tab<-matrix(0,nbprod,lastvar-firstvar+1)

if (method =="mean"){
  for (j in firstvar:lastvar){
    for (i in 1:nbprod){             
      tab[i,j-firstvar+1]<-mean(donnee[donnee[,col.p]==levels(donnee[,col.p])[i],j],na.rm=TRUE)
    }
  }
}

if (method =="coeff"){

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
    aux <- summary.lm(aov( formule, data = donnee, subset=subset,na.action =na.exclude))$coef  
    tab[-nbprod,varendo-firstvar+1] <- aux[2:nbprod,1]
    tab[nbprod,varendo-firstvar+1] <-  - sum(tab[-nbprod,varendo-firstvar+1])
    tab[,varendo-firstvar+1] <-  tab[,varendo-firstvar+1]+aux[1,1]
  }
}

dimnames(tab) = list(levels(donnee[,col.p]),labels(donnee)[[2]][firstvar:lastvar])
tab=as.data.frame(tab)
if (length(file)!=0) write.csv2(tab,file=file,sep=";",dec=",")
return(tab)

}
