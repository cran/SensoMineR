"gpa" <-
function(df,tolerance=10^-10,nbiteration=200,scal=TRUE,coord=c(1,2)) {
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

#-------------------------------------------------------------------------------
#   similarité
#-------------------------------------------------------------------------------
similarite<-function(X,Y)
{

if(dim(X)[[1]]!=dim(Y)[[1]]) stop("Les matrices ne sont pas de même dimensions")


# ON CENTER LES MATRICES
#####################################################
 Y<-scale(Y,scale=FALSE)
 X<-scale(X,scale=FALSE)
 
 y<-Y%*%procrustesbis(Y,X)$H
similari<-sum(diag(t(X)%*%y))/(sum(diag(t(X)%*%X))*sum(diag(t(y)%*%y)))^0.5

return(similari)
 }
 
 #-------------------------------------------------------------------------------
#   RV standardisé
#-------------------------------------------------------------------------------
 coeffRVs<-function(X,Y)
{

if(dim(X)[[1]]!=dim(Y)[[1]]) stop("Les matrices ne sont pas de même dimensions")
 n<-dim(X)[[1]]

# ON CENTER LES MATRICES

 Y<-scale(Y,scale=FALSE)
 X<-scale(X,scale=FALSE)
   rv<-coeffRV(X,Y)
 
# on fait un calcul exact dans le cas de 2 lignes et 3 lignes
if(n<4)
{
  perm24<-permute(n)

listX<-array(0,c(dim(X)[[1]],dim(X)[[2]],dim(perm24)[[1]]))
                            
for (i in 1:dim(perm24)[[1]])
{
listX[,,i]<- scale(X[perm24[i,],],scale=FALSE)

}

listRV<-NULL
 listRVbis<-NULL

for (i in 1:dim(perm24)[[1]])
{
listRV<-c(listRV,coeffRV(listX[,,i],Y))
#listRVbis<-c(listRVbis,sum(diag(listX[,,i]%*%t(listX[,,i])%*%Y%*%t(Y))))
}

esperance<-mean(listRV)
variance<-sum((listRV-esperance)^2)/dim(perm24)[[1]]
    

}
else
{

 
#ON CALCULES LES MATRICES DE PRODUIT SCALAIRE

betax<-(sum(diag(X%*%t(X))))^2/sum(diag(X%*%t(X)%*%X%*%t(X)))
betay<-(sum(diag(Y%*%t(Y))))^2/sum(diag(Y%*%t(Y)%*%Y%*%t(Y)))

alphax<-n-1-betax
alphay<-n-1-betay

deltax<-sum(diag(X%*%t(X))^2)/sum(diag(X%*%t(X)%*%X%*%t(X)))
gammax<-(n-1)*(n*(n+1)*deltax-(n-1)*(betax+2))/((n-3)*(n-1-betax))
deltay<-sum(diag(Y%*%t(Y))^2)/sum(diag(Y%*%t(Y)%*%Y%*%t(Y)))
gammay<-(n-1)/((n-3)*(n-1-betay))*(n*(n+1)*deltay-(n-1)*(betay+2))


esperance<-(betax^0.5)*betay^0.5/(n-1)
variance<-2*alphay*alphax/((n+1)*(n-1)^2*(n-2))*(1+(n-3)*gammax*gammay/(2*n*(n-1)))
 }
rvstd<-(rv-esperance)/variance^0.5
return(list(rvstd=rvstd,rv=rv,moyenne=esperance,variance=variance))
}


#CALCUL DES PERMUTATIONS
multiply.vec <- function(vec) {
        k <- length(vec)
        rep(vec,rep(factorial(k-1),k))}

permute <- function(nbp) {
if (nbp==2)
{ perm<-rbind(c(1,2),c(2,1))
}
else
{
perm <- matrix(nrow=factorial(nbp),ncol=nbp)
perm[,1] <- multiply.vec(1:nbp)
for (j in 2:(nbp-1)) {
restants <- apply(matrix(perm[seq(1,factorial(nbp),factorial(nbp-j+1)),1:(j-1)],
ncol=j-1,nrow=length(seq(1,factorial(nbp),factorial(nbp-j+1)))),1,setdiff,x=1:nbp)
perm[,j] <- as.vector(apply(restants,2,multiply.vec)) }
perm[,nbp] <- apply(perm[,1:nbp-1],1,setdiff,x=1:nbp)
}
perm }
#--------------------------------------------------------
# Coeff RV
#-------------------------------------------------------
coeffRV<-function(X,Y)
{

if(dim(X)[[1]]!=dim(Y)[[1]]) stop("Les matrices ne sont pas de même dimensions")


# ON CENTER LES MATRICES

 Y<-scale(Y,scale=FALSE)
 X<-scale(X,scale=FALSE)
 
#ON CALCULES LES MATRICES DE PRODUIT SCALAIRE

W1<-X%*%t(X)
W2<-Y%*%t(Y)

rv<-sum(diag(W1%*%W2))/(sum(diag(W1%*%W1))*sum(diag(W2%*%W2)))^0.5
return(rv)
}

################################################################################
# Programme permettant les sorties de PANOVA pour les tableaux contenants 
# des valeurs manquantes
################################################################################
crit.procGPAcvmqte<-function(x)
{
#x est un objet de la classe procGPAc
if (!inherits(x, "GPAc")) 
        stop("Object of type 'GPAc' expected")
contribindiv<-matrix(0,dim(x$Xfin)[[1]],3)
contriconfig<-matrix(0,dim(x$Xfin)[[3]],3)
contridim<-matrix(0,dim(x$consensus)[[2]],3)
nomligne<-row.names(x$depart)
nomconfig<-tab.names(x$depart)

Maii<-0*(x$M[,,1]%*%x$consensus%*%t(x$M[,,1]%*%x$consensus))
for (i in 1:dim(x$M)[[3]])
{
Maii<-Maii+x$M[,,i]%*%x$consensus%*%t(x$M[,,i]%*%x$consensus)


}

aii<-diag(Maii)

Meii<-x$Xfin[,,1]-x$consensus
Mbii<-0*(x$Cj[,,1]%*%Meii%*%t(x$Cj[,,1]%*%Meii))

for (i in 1:dim(x$M)[[3]])
{
Meii<-x$Xfin[,,i]-x$consensus
Mbii<-Mbii+(x$Cj[,,i]%*%Meii%*%t(x$Cj[,,i]%*%Meii))

}
bii<-diag(Mbii)
Mdii<-0*x$Xfin[,,1]%*%t(x$Xfin[,,1])
for (i in 1:dim(x$M)[[3]])
{
Mdii<-Mdii+x$Xfin[,,i]%*%t(x$Xfin[,,i])
}

dii<-diag(Mdii)


for (i in 1:dim(x$Xfin)[[1]])
{
contribindiv[i,1]<-aii[i]

contribindiv[i,2]<-bii[i]
contribindiv[i,3]<-dii[i]
}

contribindiv<-rbind(contribindiv,colSums(contribindiv))
rownames(contribindiv)<-c(nomligne,"sum")
colnames(contribindiv)<-c("SSfit","SSresidual","SStotal")

contriconfig<-matrix(0,dim(x$Xfin)[[3]],3)
contridim<-matrix(0,dim(x$consensus)[[2]],3)
for (i in 1:dim(x$M)[[3]])
{



Meii<-x$Xfin[,,i]-x$consensus
contriconfig[i,1]<-x$poids[i]*sum(diag(t(x$Z)%*%x$Cj[,,i]%*%x$Xdeb[,,i]%*%x$R[,,i]))
contriconfig[i,2]<-sum(diag(t(Meii)%*%x$Cj[,,i]%*%Meii))
contriconfig[i,3]<-x$poids[i]^2*(sum(diag(t(x$Xdeb[,,i])%*%x$Cj[,,i]%*%x$Xdeb[,,i])))-sum(diag(t(x$consensus)%*%x$Cj[,,i]%*%Meii))

}
contriconfig<-rbind(contriconfig,colSums(contriconfig))
rownames(contriconfig)<-c(nomconfig,"sum")
colnames(contriconfig)<-c("SSfit","SSresidual","SStotal")


###########################################################
# Décomposition par dimension p49
#############################################################
Mfii<-0*(t(x$Xfin[,,1]-x$consensus)%*%x$Cj[,,1]%*%(x$Xfin[,,1]-x$consensus))
Mgii<-0*(t(x$poids[1]*x$Cj[,,1]%*%x$Xdeb[,,1]%*%x$R[,,1]%*%x$K)%*%(x$poids[1]*x$Cj[,,1]%*%x$Xdeb[,,1]%*%x$R[,,1]%*%x$K))
for (i in 1:dim(x$M)[[3]])
{
Mfii<-Mfii+t(x$Xfin[,,i]-x$consensus)%*%x$Cj[,,i]%*%(x$Xfin[,,i]-x$consensus)
Mgii<-Mgii+t(x$poids[i]*x$Cj[,,i]%*%x$Xdeb[,,i]%*%x$R[,,i]%*%x$K)%*%(x$poids[i]*x$Cj[,,i]%*%x$Xdeb[,,i]%*%x$R[,,i]%*%x$K)
}

for (i in 1:dim(x$consensus)[[2]])
{

contridim[i,1]<-diag(x$gama)[i]
contridim[i,2]<-diag(Mfii)[i]
contridim[i,3]<-diag(Mgii)[i]

}
contridim<-rbind(contridim,colSums(contridim))

rownames(contridim)<-c(c(1:dim(x$consensus)[[2]]),"sum")
colnames(contridim)<-c("SSfit","SSresidual","SStotal")



contribution<-list()
contribution$objet<-contribindiv
contribution$config<-contriconfig
contribution$dim<-contridim

return(contribution)
}

################################################################################
# Programme permettant les sorties de PANOVA pour les tableaux sans valeurs 
#  manquantes
################################################################################

crit.procGPAcsansvm<-function(x)
{
  #x est un objet de la classe procGPAc
if (!inherits(x, "GPAc")) 
        stop("Object of type 'GPAc' expected")
      
contribindiv<-matrix(0,dim(x$Xfin)[[1]],3)
contribconfig<-matrix(0,dim(x$Xfin)[[3]],3)
contridim<-matrix(0,dim(x$consensus)[[2]],3)
nbj<-dim(x$Xfin)[[3]]
nomligne<-row.names(x$depart)
nomconfig<-tab.names(x$depart)
s<-NULL
s1<-NULL
s3<-NULL
for(i in 1:dim(x$Xfin)[[3]])
        {
            
            s<-rbind(s,colSums(x$Xfin[,,i]^2,na.rm=FALSE,dims=1))
            s1<-rbind(s1,colSums((x$Xfin[,,i]-x$consensus)^2,na.rm=FALSE,dims=1))
             s3<-rbind(s1,colSums((x$consensus)^2,na.rm=FALSE,dims=1))
        }
      
        contribconfig[,3]<-rowSums(s,na.rm=FALSE,dims=1)
        contribconfig[,2]<-rowSums(s1,na.rm=FALSE,dims=1)
        contribconfig[,1]<-0
       
        contribconfig<-rbind(contribconfig,colSums(contribconfig))

        colnames(contribconfig)<-c("SSfit","SSresidual","SStotal")
        row.names(contribconfig)<-c(nomconfig,"sum")
     
        
        
        
        s<-0*rowSums(x$Xfin[,,1]^2,na.rm=FALSE,dims=1)
        s1<-0*rowSums((x$Xfin[,,i]-x$consensus)^2,na.rm=FALSE,dims=1)
        
for(i in 1:dim(x$Xfin)[[3]])
        {
            
            s<-s+rowSums(x$Xfin[,,i]^2,na.rm=FALSE,dims=1)
            s1<-s1 +rowSums((x$Xfin[,,i]-x$consensus)^2,na.rm=FALSE,dims=1)
             
        }
        
        contribindiv[,3]<-s
        contribindiv[,2]<-s1
        contribindiv[,1]<-nbj*rowSums((x$consensus)^2,na.rm=FALSE,dims=1)
       
        contribindiv<-rbind(contribindiv,colSums(contribindiv))

        colnames(contribindiv)<-c("SSfit","SSresidual","SStotal")
        row.names(contribindiv)<-c(nomligne,"sum")
 s<-NULL
 s1<-NULL
 s3<-NULL
  for(i in 1:dim(x$Xfin)[[3]])
        {
            
            s<-rbind(s,colSums(x$Xfin[,,i]^2,na.rm=FALSE))
            s1<-rbind(s1,colSums((x$Xfin[,,i]-x$consensus)^2,na.rm=FALSE,dims=1))
             s3<-rbind(s1,colSums((x$consensus)^2,na.rm=FALSE,dims=1))
        }   
                     s3<-colSums((x$consensus)^2,na.rm=FALSE,dims=1)
                
                
contridim<-cbind( nbj*(s3),(colSums(s1)),(colSums(s)))
contridim<-rbind(contridim,colSums(contridim))

  colnames(contridim)<-c("Consensus","residus","Total")
        row.names(contridim)<-c(paste("dim",1:dim(x$Xfin)[[2]]),"Total")
       
contribution<-list()
contribution$objet<-contribindiv/nbj*100
contribution$config<-contribconfig/nbj*100
contribution$dimension<- contridim /nbj*100
return(contribution)

}
########################################################
"summaryGPAc" <-
function(x) {
    if (!inherits(x, "GPAc")) stop("Object of type 'GPAc' expected")
    lab <- c("M","Cj","Xdeb","it","poids","translation","Z","K","gama","Xfin","consensus","R","VMQTE")
    tab <- c(paste(dim(x$M)[1],"x",dim(x$M)[2],"x",dim(x$M)[3]),"array", " positions NA")
    tab <- rbind(tab,c( paste(dim(x$Cj)[[1]],"x",dim(x$Cj)[[2]],"x",dim(x$Cj)[[3]]),"array", "aide centrage") )
    tab <- rbind(tab,c(paste(dim(x$Xdeb)[1],"x",dim(x$Xdeb)[2],"x",dim(x$Xdeb)[3]),"array", "matrices init. NA recod.") )
    tab <- rbind(tab,c( paste(dim(x$it)[1],"x",dim(x$it)[2]),"matrice", " Crit. après chaque itération"))
    tab <- rbind(tab,c( paste(dim(x$poids)[[1]],"x",dim(x$poids)[[2]]),"vecteur", "Vecteur des poids"))
    tab <- rbind(tab,c( paste(dim(x$translation)[[1]],"x",dim(x$translation)[[2]]),"matrice", "Vect. trans. pr chaque ss config."))
    tab <- rbind(tab,c(paste(dim(x$Z)[1],"x",dim(x$Z)[[2]]),"matrice", "Z avant proj."))
    tab <- rbind(tab,c( paste(dim(x$K)[1],"x",dim(x$K)[2]),"matrice", "mat de proj."))
    tab <- rbind(tab,c( paste(dim(x$gama)[1],"x",dim(x$gama)[2]),"matrice", "mat des vp"))
    tab <- rbind(tab,c( paste(dim(x$Xfin)[1],"x",dim(x$Xfin)[2],"x",dim(x$Xfin)[3]),"array", "Config. ind. fin."))
    tab <- rbind(tab,c( paste(dim(x$consensus)[1],"x",dim(x$consensus)[2]),"matrice", "Config. moyenne fin."))
    tab <- rbind(tab,c( paste(dim(x$R)[1],"x",dim(x$R)[2],"x",dim(x$R)[3]),"array", "Transf. othogonales"))
    tab <- rbind(tab,c(length(x$VMQTE),"Booleen", "Présence ou pas de NA"))
    row.names(tab) <- lab
    tab <- as.data.frame(tab)
    tab
}
########################################################
    "placevm" <- function(mat) {
    if(!is.matrix(mat))  stop("A matrix please !!")
    vect <- NULL
    matri <- NULL
    i <- 1
    while (i<=dim(mat)[[1]]) {
#if(sum(is.na(mat[i,]))==dim(mat)[[2]]) {
if (any(is.na(mat[i,]))) {
    matri <- matri
            vect <- c(vect,i)
            i <- i+1
        }
    else {
    matri <- rbind(matri,mat[i,])
    i <- i+1
        }
    }
    sol <- list()
    sol$matri <- matri
    sol$vect <- vect
    return(sol)
}
########################################################
"procrustesbis" <-
function(X1,X2) {
    if(!is.matrix(X1))  stop("On souhaite ici avoir un tableau de type matrice  ")
    if(!is.matrix(X2))  stop("On souhaite ici avoir un tableau de type matrice  ")
    if(dim(X2)[[2]]!=dim(X1)[[2]]) stop("On souhaite ici avoir un tableau de type matrice ayant le même nombre de colonnes  ")
    nbcolonne <- dim(X1)[[2]]
    X1c <- scale(X1,scale=FALSE)
    X2c <- scale(X2,scale=FALSE)
    Aj <- t(X1c)%*%X2c
    lola <- eigen(t(Aj)%*%Aj)
    Qj <- as.matrix(eigen(t(Aj)%*%Aj)$vectors)
    phij <- diag (lola$values,length(lola$values),length(lola$values))
    if (sum(lola$values>10^-15)==length(lola$values)) {
 Pj <- Aj%*%Qj%*%as.matrix(diag(diag(phij)^(-0.5),length(lola$values),length(lola$values)))
        H <- Pj%*%t(Qj)
    }
    else {
        nbracinepos <- as.numeric(sum((lola$values)>10^-15))
        Pjstar <- Aj%*%Qj[,1:nbracinepos]%*%as.matrix(diag((diag(phij)[1:nbracinepos])^(-0.5),length(lola$values[1:nbracinepos]),length(lola$values[1:nbracinepos])))
        Pbar <- matrix(0,nbcolonne,(length(lola$values)-nbracinepos))
        Pjstar1 <- Pjstar
        for (k in 1:(length(lola$values)-nbracinepos)) {
            yinit <- rnorm(nbcolonne)
            yinit.lm <- lm(yinit~0+Pjstar1)
            vectorth <- yinit-fitted.values(yinit.lm)
            Pbar[,k] <- vectorth/(t(vectorth)%*%vectorth)^0.5
            if (sum (sign(Pbar[,k]/Qj[,(nbracinepos+k)]),na.rm=TRUE)<0) {
           Pjstar1 <- cbind(Pjstar1,-Pbar[,k])
            }
            else {
            Pjstar1 <- cbind(Pjstar1,Pbar[,k])
            }
        }
        H <- Pjstar1%*%t(as.matrix(lola$vectors))
    }
    rho <- sum(diag(X1c%*%H%*%t(X2c)))/sum(diag(X1c%*%t(X1)))
    result <- list()
    result$H <- H
    result$rho <- rho
    return(result)
}
########################################################

    if(!is.ktab(df)) stop("not ktab")
    blo <- df$blo
    nbj <- length(blo)
    taille <- max(blo)
    X <- array(0,c(dim(df[[1]])[[1]],taille,nbj))
    for (i in 1:nbj) {
X[,1:blo[i],i] <- as.matrix(df[[i]])
    }
    Xm <- NULL
    v <- NULL
    M <- NULL
    C <- NULL
    U <- NULL
    Ip <- NULL
    vdiag <- NULL
    mat1 <- NULL
    Cc <- NULL
    invgC <- NULL
    Xm1 <- NULL
    nbj <- NULL
    nbjuge <- NULL
    pds <- NULL
    nbjuge <- dim(X)[[3]]
    p <- dim(X)[[1]]
    nbcolonne <- dim(X)[[2]]
    M <- array(0,c(p,p,nbjuge))
    Cj <- array(0,c(p,p,nbjuge))
    U <- matrix(1,p,1)    
    Ip <- diag(rep(1,p),p,p)
    Xm <- X
    v <- rep(1,p)
    VMQTE <- FALSE
    for (j in 1:nbjuge) {
vdiag <- v
    if (length(placevm(Xm[,,j])$vect)==0) {
        M[,,j] <- Ip
        }
        else {
        vdiag[placevm(Xm[,,j])$vect] <- 0
        M[,,j] <- diag(vdiag,length(vdiag),length(vdiag))
        VMQTE <- TRUE
        }
    }
    vdiag <- NULL
    mat1 <- U%*%t(U)
    for (j in 1:nbjuge) {
    Cj[,,j] <- M[,,j]%*%(Ip-mat1%*%M[,,j]/as.numeric(t(U)%*%M[,,j]%*%U))
    }
    Cc <- Cj[,,1]
    for (j in 2:nbjuge) {
    Cc <- Cc+Cj[,,j]
    }
    invgC <- ginv(Cc)
    Xm1 <- Xm
    for (j in 1:nbjuge) {
    for (i in 1:dim(Xm)[[1]]) {
        Xm1[i,,j] <- replace(Xm1[i,,j],is.na(Xm1[i,,j]),999999)
    }
    }
    lambda2 <- 0
    for (j in 1:nbjuge) {
    lambda2 <- lambda2+sum(diag(t(Xm1[,,j])%*%Cj[,,j]%*%Xm1[,,j]))
    }
    lambda <- (nbjuge/(lambda2))^0.5
    Xnorm <- Xm1*lambda
    diagW <- NULL
    for (i in 1:nbjuge) {
    diagW <- c(diagW,1/sum(diag(t(Xnorm[,,i])%*%Cj[,,i]%*%Xnorm[,,i]))^0.5)
    }
    W12 <- diag(diagW,nbjuge,nbjuge)
    pds <- rep(1,nbjuge)
    R <- array(0,c(nbcolonne,nbcolonne,nbjuge))
    Im <- diag(rep(1,nbcolonne),nbcolonne,nbcolonne)
    for (i in 1:nbjuge) {
    R[,,i] <- Im 
    }
    Aj <- array(0,c(nbcolonne,nbcolonne,nbjuge))
    for (j in 1 : nbjuge) {   
sommetemp <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
    for (i in 2:nbjuge) {
        sommetemp <- sommetemp+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
    }
    R[,,j] <- procrustesbis(Cj[,,j]%*%Xnorm[,,j],invgC%*%(sommetemp-pds[j]*Cj[,,j]%*%Xnorm[,,j]%*%R[,,j]))$H
    }
    sommetemp2 <- NULL
    sommetemp2 <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
    for (i in 2:nbjuge) {
    sommetemp2 <- sommetemp2+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
    }
    matidd <- t(sommetemp2)%*%invgC%*%sommetemp2
    lossf <- nbjuge-sum(diag(t(sommetemp2)%*%invgC%*%sommetemp2))
    lossf2 <- lossf
    if(scal) {
matY <- matrix(0,nbjuge,nbjuge)
B <- array(0,c(p,nbcolonne,nbjuge))
    for (k in 1:nbjuge) {
    B[,,k] <- Cj[,,k]%*%Xnorm[,,k]%*%R[,,k]
        }
for (k in 1:nbjuge) {
        for ( l in 1 :nbjuge) {
matY[k,l] <- sum(diag(t(B[,,k])%*%invgC%*%B[,,l]))
        }
    }
eigzou <- eigen(W12%*%matY%*%W12)
verifsigne <- sum(eigzou$vectors[,1]<0)
tailleeig <- dim(eigzou$vectors)[[1]]
if(verifsigne==tailleeig) {
    vecteurpropre <- eigzou$vectors[,1]*(-1)
}
else {
    vecteurpropre<-eigzou$vectors[,1]
}
verifsigne <- NULL
tailleeig <- NULL
pds <- (nbjuge)^0.5*W12%*%as.matrix(vecteurpropre)
sommetemp2 <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
for (i in 2:nbjuge) {
    sommetemp2 <- sommetemp2+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
    }
lossf2 <- (nbjuge-sum(diag(t(sommetemp2)%*%invgC%*%sommetemp2)))
    }
    tol <- lossf2
    itorth <- lossf
    itpoids <- lossf2
    lossf <- lossf2
    compteur <- 0
    while (tol>tolerance && compteur<nbiteration) {
    for (j in 1 : nbjuge) {   
            sommetemp <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
    for (i in 2:nbjuge) {
            sommetemp <- sommetemp+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
            }
            R[,,j] <- procrustesbis(Cj[,,j]%*%Xnorm[,,j],invgC%*%(sommetemp-pds[j]*Cj[,,j]%*%Xnorm[,,j]%*%R[,,j]))$H
        }
sommetemp2 <- NULL
sommetemp2 <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
for (i in 2:nbjuge) {
        sommetemp2 <- sommetemp2+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
    }
matidd <- t(sommetemp2)%*%invgC%*%sommetemp2
lossf2ortho <- nbjuge-sum(diag(t(sommetemp2)%*%invgC%*%sommetemp2))
lossf2 <- lossf2ortho
lossfpoids <- 0
    if (scal) {
    matY <- matrix(0,nbjuge,nbjuge)
        B <- array(0,c(p,nbcolonne,nbjuge))
        for (k in 1:nbjuge) {
        B[,,k] <- Cj[,,k]%*%Xnorm[,,k]%*%R[,,k]
            }
    for (k in 1:nbjuge) {
    for ( l in 1 :nbjuge) {
        matY[k,l] <- sum(diag(t(B[,,k])%*%invgC%*%B[,,l]))
    }
        }
    eigzou <- eigen(W12%*%matY%*%W12)
    if(sum(eigzou$vectors[,1]<0)==dim(eigzou$vectors)[[1]]) {
vecteurpropre <- -eigzou$vectors[,1]
    }
    else {
vecteurpropre <- eigzou$vectors[,1]
    }
    pds <- (nbjuge)^0.5*W12%*%as.matrix(vecteurpropre)
    sommetemp2 <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
    for (i in 2:nbjuge) {
        sommetemp2 <- sommetemp2+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
        }
    lossfpoids <- (nbjuge-sum(diag(t(sommetemp2)%*%invgC%*%sommetemp2)))
    lossf2 <- (nbjuge-sum(diag(t(sommetemp2)%*%invgC%*%sommetemp2)))
}
tol <- (lossf-lossf2)
lossf <- lossf2
itorth <- c(itorth,lossf2ortho)
itpoids <- c(itpoids,lossfpoids)
compteur <- compteur+1
    }
    sommetemp2 <- NULL
    sommetemp2 <- pds[1]*Cj[,,1]%*%Xnorm[,,1]%*%R[,,1]
    for (i in 2:nbjuge) {
    sommetemp2 <- sommetemp2+pds[i]*Cj[,,i]%*%Xnorm[,,i]%*%R[,,i]
    }
    pp <- invgC%*%sommetemp2
    translation <- matrix(0,nbcolonne,nbjuge)
    for (i in 1:nbjuge) {
translation[,i] <- (t(pds[i]*Xnorm[,,i]-pp%*%t(R[,,i]))%*%M[,,i]%*%U)/as.numeric(pds[i]*t(U)%*%M[,,i]%*%U)
    }
    ppeig <- eigen(t(pp)%*%Cc%*%pp)
    Xfin <- Xnorm
    for (k in 1 :nbjuge) {
Xfin[,,k] <- pds[k]*M[,,k]%*%(Xnorm[,,k]-U%*%t(translation[,k]))%*%R[,,k]%*%as.matrix(ppeig$vectors)
    }
    it <- cbind(itorth,itpoids)
    colnames(it) <- c("rotation step","scaling step")
        par()
        res <- Xfin[, , 1]
        for (i in 2:dim(Xfin)[3]) {
            res <- rbind(res, Xfin[, , i])
        }
        coo.part <- data.frame(res[, coord])
        row.names(coo.part) <- paste(rep(1:dim(Xfin)[1],dim(Xfin)[3]), rep(1:dim(Xfin)[3],
            each = dim(Xfin)[1]), sep = ".")
        fact.gpa <- factor(rep((1:dim(Xfin)[1]), dim(Xfin)[3]))
        consensus <- data.frame(pp%*%as.matrix(ppeig$vectors))
        lab.gpa <- row.names(consensus)
        senso.class(coo.part, fac = fact.gpa)
        title(sub = paste("Row projection (", "comp", coord[1],
            "-", "comp", coord[2], ")"), adj = 0.5, line = 3.8)

      consensus<-pp%*%as.matrix(ppeig$vectors)
# calcul de la dimension max du consensus
cte<-1

while(cte<=dim(consensus)[[2]])
{
if(all(abs(consensus[,cte])<10^-15))
{
fin<-cte-1
cte=dim(consensus)[[2]]+1
}
else
{
fin<-cte
cte<-cte+1

}

}



      row.names(consensus)<-row.names(df)
      row.names(Xfin)<-row.names(df)
    result <- list()
    class(result) <- c("GPAc", "list")
    result$call <- match.call()
     result$depart<-df
    result$M <- M
    result$Cj <- Cj
    result$consensus <-consensus[,1:fin]
    result$Z <- pp
    result$Xfin <- Xfin[,1:fin,]
    result$Xdeb <- Xnorm
    result$poids <- pds
    result$translation <- translation
    result$it <- it
    result$Rj <- R
    result$K <- as.matrix(ppeig$vectors)
    result$gama <- as.matrix(diag(ppeig$values))
    result$VMQTE <- VMQTE
    result$resume <- summaryGPAc(result)
   
  
##############################################################
nbjuge<-length(df$blo)
x<-result

# on calcule les coefficients RV
     RVs<-matrix(-1,nbjuge,nbjuge)
        RV<-matrix(-1,nbjuge,nbjuge)
        sim<-matrix(-1,nbjuge,nbjuge)           


if(x$VMQTE)
{
x<-x
X<-array(0,c(dim(df[[1]])[[1]],max(df$blo),nbjuge))
for (i in 1:nbjuge)
{
X[,1:df$blo[[i]],i]<-as.matrix(df[[i]])
}


# on récupère les emplacements des valeurs manquantes
vmplacelist<-list()
length(vmplacelist)<-nbjuge

for (theta in 1:nbjuge)
{
if (length(placevm(X[,,theta])$vect)!=0)
{
vmplacelist[[theta]]<-placevm(X[,,theta])$vect
}


}

for (i in 1:nbjuge)
{
for (j in 1:nbjuge)
{
if (length(c(vmplacelist[[i]],vmplacelist[[j]]))!=0){
Xi<-X[,,i][-c(vmplacelist[[i]],vmplacelist[[j]]),]
Xj<-X[,,j][-c(vmplacelist[[i]],vmplacelist[[j]]),]
}

if (length(c(vmplacelist[[i]],vmplacelist[[j]]))==0){
Xi<-X[,,i]
Xj<-X[,,j]
}

          RVs[i,j]<-coeffRVs(Xi,Xj)$rvstd
         RV[i,j]<-coeffRVs(Xi,Xj)$rv
         sim[i,j]<-similarite(Xi,Xj)
        
        
}
}

}
else
{

X<-array(0,c(dim(df[[1]])[[1]],max(df$blo),nbjuge))
for (i in 1:nbjuge)
{
X[,1:df$blo[[i]],i]<-as.matrix(df[[i]])
}


for (i in 1:nbjuge)
{
for (j in 1:nbjuge)
{
Xi<-X[,,i]
Xj<-X[,,j]


        RVs[i,j]<-coeffRVs(Xi,Xj)$rvstd
         RV[i,j]<-coeffRVs(Xi,Xj)$rv
         sim[i,j]<-similarite(Xi,Xj)
        
}
}  
}
 
##############################################################
   row.names(RV)<-colnames(RV)<-tab.names(df)
   row.names(RVs)<-colnames(RVs)<-tab.names(df)
   row.names(sim)<-colnames(sim)<-tab.names(df)
    resultat<-list()
    resultat$RV<-RV
    resultat$RVs<-RVs
    resultat$simi<-sim
    resultat$scaling<-result$poids
    resultat$consensus<-result$consensus
    resultat$Xfin<-result$Xfin

    if (result$VMQTE) resultat$PANOVA<-crit.procGPAcvmqte(result)
    else resultat$PANOVA<-crit.procGPAcsansvm(result)
    
    return(resultat)
}
