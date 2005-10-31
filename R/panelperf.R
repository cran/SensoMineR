"panelperf" <- function(donnee,formul,subset=NULL,firstvar,lastvar=ncol(donnee),random=TRUE){

  for (j in 1 :(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])
  formul = as.formula(formul)
  lab <- labels(donnee)[[2]]
  equation <- as.character(formul)

  dim.donnee <- ncol(donnee)
  for (i in 1:dim.donnee) {
    if (gsub(" ","",strsplit(equation,split="+",extended=FALSE)[[2]][1])==lab[i]) col.p <- i
    if (gsub(" ","",strsplit(equation,split="+",extended=FALSE)[[2]][2])==lab[i]) col.j <- i
  }
  res <- matrix(0,lastvar+1-firstvar,1)
  r2 <- matrix(0,lastvar+1-firstvar,1)
  variab <- perf <- matrix(0,lastvar+1-firstvar,length(strsplit(equation[2],"+",extended=FALSE)[[1]]))

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
    aux1 <- aov( formule , data = donnee, subset=subset,na.action =na.exclude)
    aux <- summary(aux1)[[1]]
    perf[varendo-firstvar+1,] <- aux[-nrow(aux),5]
    variab[varendo-firstvar+1,] <- aux[-nrow(aux),2]/sum(aux[,2])
    res[varendo-firstvar+1,] <- sqrt(aux[nrow(aux),3])
    r2[varendo-firstvar+1,] <- summary.lm(aux1)$r.squared
    
    if (random) {
      for (i in 1:nrow(aux)){
       if (gsub(" ","",paste("C(",colnames(donnee)[col.p],",sum):C(",colnames(donnee)[col.j],",sum)")) == gsub(" ","",rownames(aux)[i])) row.interact=i
       if (gsub(" ","",paste("C(",colnames(donnee)[col.j],",sum):C(",colnames(donnee)[col.p],",sum)")) == gsub(" ","",rownames(aux)[i])) row.interact=i
      }
      perf[varendo-firstvar+1,1] <- pf(aux[1,3]/aux[row.interact,3],aux[1,1],aux[row.interact,1],lower.tail=FALSE)
    }
  }

aa <- strsplit(as.character(formule),split="~",extended=FALSE)[[3]]
aa <- gsub("C\\(","",aa)
aa <- gsub(", sum","",aa)
aa <- gsub(")","",aa)
dimnames(perf) <- list(labels(donnee)[[2]][firstvar:lastvar],gsub(" ","",gsub(", sum)","",gsub("C\\(","",rownames(aux))))[-nrow(aux)])
dimnames(variab) <- list(labels(donnee)[[2]][firstvar:lastvar],gsub(" ","",gsub(", sum)","",gsub("C\\(","",rownames(aux))))[-nrow(aux)])
dimnames(res) <- list(labels(donnee)[[2]][firstvar:lastvar],"stdev residual")
dimnames(r2) <- list(labels(donnee)[[2]][firstvar:lastvar],"r2")

panelperf = list() 
panelperf$p.value = perf
panelperf$variability = variab
panelperf$res = res
panelperf$r2 = r2
return(panelperf)
}
