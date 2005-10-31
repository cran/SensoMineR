"coltable" <-function(matrice,nbrow=nrow(matrice),nbcol=ncol(matrice),level.lower=0.05,col.lower="mistyrose",level.upper=1e10,col.upper="lightblue",cex=0,nbdec=4,main.title=NULL) {

################################################################
"fill" <- function(matrice,nbrow,nbcol,pol,level.lower,col.lower="mistyrose",level.upper,col.upper="lightblue",main.title=NULL){

#cadre
dim1 <- dim(matrice)[1] 
dim2 <- dim(matrice)[2] 
for (i in 0:dim1) rect(0,1-i*(1/(nbrow+1)),1/(nbcol+1),1-(i+1)*(1/(nbrow+1)),col="white",border=NULL)
for (j in 1:(dim2-1)) rect(j*(1/(nbcol+1)),1,(j+1)*(1/(nbcol+1)),1-(1/(nbrow+1)),col="white",border=NULL)

for (j in 1:(dim2-1)){ for (i in 1:dim1){

if (is.na(matrice[i,j+1])){
   rect(j*(1/(nbcol+1)),1-i*(1/(nbrow+1)),(j+1)*(1/(nbcol+1)),1-(i+1)*(1/(nbrow+1)),col="gray",border=NULL)
}

else { if (matrice[i,j+1]<level.lower){
         rect(j*(1/(nbcol+1)),1-i*(1/(nbrow+1)),(j+1)*(1/(nbcol+1)),1-(i+1)*(1/(nbrow+1)),col=col.lower,border=NULL)
       }
       else { if (matrice[i,j+1]>level.upper){
            rect(j*(1/(nbcol+1)),1-i*(1/(nbrow+1)),(j+1)*(1/(nbcol+1)),1-(i+1)*(1/(nbrow+1)),col=col.upper,border=NULL)
            }
            else {
             rect(j*(1/(nbcol+1)),1-i*(1/(nbrow+1)),(j+1)*(1/(nbcol+1)),1-(i+1)*(1/(nbrow+1)),col="white",border=NULL)
           }
       }    
   }
 }
}

#fill
    dim1 <- dim(matrice)[1]
    dim2 <- dim(matrice)[2]
    for (i in 1:dim1){
    for (j in 0:(dim2-1)) {
        text((j+0.5)*(1/(nbcol+1)),1-(i+0.5)*(1/(nbrow+1)),matrice[i,j+1],cex=pol)
    }
    }
#titre

    for (j in 0:nbcol) {
    text((j+0.5)*(1/(nbcol+1)),1-(1/(nbrow+1))/2,names(matrice)[j+1],cex=pol)
    }

}
################################################################

################################################################
"police" <- function(matrice,nbrow,nbcol,nbdec) {
    get(getOption("device"))(12,8)
    par(mar=c(0,0,2,0))
    plot.new() ; title(main=main.title);
    a <- c(rownames(matrice),colnames(matrice))
    nb=NULL
    for (i in 1:nbdec) nb <- paste(nb,"0",sep="")
    nb <- paste(nb,"0.e-00")
    a <- c(a,nb)
    b <- min(nbcol,15)
    return((round((1/(b+1))/max(strwidth(a)),2)*100-5)/100)
}
################################################################

matrice <- signif(matrice,nbdec)
matrice=cbind.data.frame(rownames(matrice),matrice)
colnames(matrice)[1]=" "
dim1 <- nrow(matrice)
dim2 <- ncol(matrice)
dim2 <- dim2-1
size <- cex
if (nbrow>dim1){ nbrow <- dim1 }
if (nbcol>dim2){ nbcol <- dim2 }
if (dim2%/%nbcol==dim2/nbcol) {
 for (j in 0:(dim2%/%nbcol-1)) {
    for (i in 0:(dim1%/%nbrow-1)){
        A <- data.frame(matrice[(i*nbrow+1):((i+1)*nbrow),1])
        names(A)=names(matrice)[1]
        B <- matrice[(i*nbrow+1):((i+1)*nbrow),(1+j*nbcol+1):(1+(j+1)*nbcol)]
        B <- cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
          get(getOption("device"))(12,8)
          par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
        fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)
    }
    if ((dim1%/%nbrow)*nbrow != dim1){
      A<-data.frame(matrice[(dim1%/%nbrow*nbrow+1):dim1,1])
      names(A)=names(matrice)[1]
      B<-data.frame(matrice[(dim1%/%nbrow*nbrow+1):dim1,(1+j*nbcol+1):(1+(j+1)*nbcol)])
      names(B)=names(matrice)[(1+j*nbcol+1):(1+(j+1)*nbcol)]
      B<-cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
         get(getOption("device"))(12,8)
         par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
      fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)

    }
  }
}
else {
    for (j in 0:(dim2%/%nbcol-1)){ #blocs de descripteurs entiers
      for (i in 0:(dim1%/%nbrow-1)){ #blocs de juges entiers
        A<-data.frame(matrice[(i*nbrow+1):((i+1)*nbrow),1])
        names(A)=names(matrice)[1]
        B<-matrice[(i*nbrow+1):((i+1)*nbrow),(1+j*nbcol+1):(1+(j+1)*nbcol)]
        B<-cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
          get(getOption("device"))(12,8)
          par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
        fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)
      }
      if ((dim1%/%nbrow)*nbrow != dim1){
        A<-data.frame(matrice[(dim1%/%nbrow*nbrow+1):dim1,1])
        names(A)=names(matrice)[1]
        B<-matrice[(dim1%/%nbrow*nbrow+1):dim1,(1+j*nbcol+1):(1+(j+1)*nbcol)]
        B<-cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
          get(getOption("device"))(12,8)
          par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
        fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)
      }
    }
    for (i in 0:(dim1%/%nbrow-1)){#pour les blocs d'individus entiers les variables qui manquent
      A<-data.frame(matrice[(i*nbrow+1):((i+1)*nbrow),1])
      names(A)=names(matrice)[1]
      B<-matrice[(i*nbrow+1):((i+1)*nbrow),(1+dim2%/%nbcol*nbcol):dim2+1]
      if (is.null(dim(B)))    B<-data.frame(B)
      names(B)=names(matrice)[(1+dim2%/%nbcol*nbcol):dim2+1]
      B<-cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
          get(getOption("device"))(12,8)
          par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
      fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)
    }
    if ((dim1%/%nbrow)*nbrow != dim1){
      A<-data.frame(matrice[(dim1%/%nbrow*nbrow+1):dim1,1]) #les individus qui manquent et les variables qui manquent
      names(A)=names(matrice)[1]
      B<-data.frame(matrice[(dim1%/%nbrow*nbrow+1):dim1,(1+dim2%/%nbcol*nbcol):dim2+1])
      if (is.null(dim(B)))    B<-data.frame(B)
      names(B)=names(matrice)[(1+dim2%/%nbcol*nbcol):dim2+1]
      B<-cbind(A,B)
        if (size==0) size <- police(matrice,nbrow,nbcol,nbdec)
        else{
          get(getOption("device"))(12,8)
          par(mar=c(0,0,2,0))
          plot.new() ; title(main=main.title);
        }
      fill(B,nbrow,nbcol,size,level.lower,col.lower,level.upper,col.upper)
    }
  }
}
