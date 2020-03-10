CA_JAR <- function(x, col.p, col.j, col.pref, jarlevel="jar"){

  ind.jar <- (1:ncol(x))[-c(col.p,col.j,col.pref)]
  nbjar <- length(ind.jar)
  for (i in ind.jar){
    x[,i]=relevel(x[,i],jarlevel)
    levels(x[,i])[1]=paste(colnames(x)[i],jarlevel,"_")
  }
  x[,col.pref]=as.factor(x[,col.pref])
  modtot <- unlist(lapply(x[,c(ind.jar,col.pref)],nlevels))
  Frequency <- matrix(NA,nrow=sum(modtot),ncol=nlevels(x[,col.p]))
  for (j in 1:ncol(Frequency)) Frequency[,j] <- unlist(lapply(x[x[,col.p]==levels(x[,col.p])[j],c(ind.jar,col.pref)],table))
  rownames(Frequency)=unlist(lapply(x[,c(ind.jar,col.pref)],levels))
  colnames(Frequency)=levels(x[,col.p])
  res <- CA(Frequency,row.sup=(sum(modtot[-length(modtot)])+1):sum(modtot), graph = FALSE)
  print(plot(res, cex=0.8, col.row=rep(3:(1+length(modtot)),modtot[-length(modtot)]),col.row.sup=rep("black",modtot[length(modtot)]),title="CA on the table product x jar variables"))
  return(list(Frequency=Frequency,res.CA=res))
}
