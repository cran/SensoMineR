"hmfa" <- function(X, H, num = 6, coord = c(1,2)) {


###############################################################"
maxmin <- function(inter,H) {

nbgroupinter <- sum(H[[2]])
nbindiv <- dim(inter)[1]/nbgroupinter

mm <- matrix(0,nbindiv,6)

for (i in 1:nbindiv) {
    mm[i,1] <- inter[i,1]
    mm[i,2] <- inter[i,1]
    mm[i,3] <- inter[i,2]
    mm[i,4] <- inter[i,2]
}

for (g in 2:nbgroupinter) {
    for (i in 1:nbindiv) {
    if (inter[(((g-1)*nbindiv)+i),1] < mm[i,1]) mm[i,1] <- inter[(((g-1)*nbindiv)+i),1] 
    if (inter[(((g-1)*nbindiv)+i),1] > mm[i,2]) mm[i,2] <- inter[(((g-1)*nbindiv)+i),1]
    if (inter[(((g-1)*nbindiv)+i),2] < mm[i,3]) mm[i,3] <- inter[(((g-1)*nbindiv)+i),2]
    if (inter[(((g-1)*nbindiv)+i),2] > mm[i,4]) mm[i,4] <- inter[(((g-1)*nbindiv)+i),2]
    }
}
return(mm)
}
###############################################################"

###############################################################"
partial.ind.tab <- function(X,H,coord=c(1,2)) {
    nbind <- dim(X)[1]
    nbnivo <- length(H)
    Xdes <- htabdes(H)
    res1 <- list()
    for (h in 1:nbnivo) {
    nbgroup <- length(H[[h]])
    res2 <- list()
    for (g in 1:nbgroup) {
        inter <- matrix(0,dim(X)[1],dim(X)[2],dimnames=list(row.names(X),variable.names(X)))
        inter1 <- ktab.data.frame(data.frame(inter),Xdes[[h]])
        inter2 <- ktab.data.frame(data.frame(X),Xdes[[h]])
        inter1[[g]] <- inter2[[g]]
        a <- inter1[[1]]
        for (g1 in 2:nbgroup) a <- cbind(a,inter1[[g1]])               
        res2[[g]] <- a
        }
    res1[[h]] <- res2
    }
    poids <- hweight(X,H)
    dilat <- hdil(H)
    interres <- dudi.pca(X,col.w=poids[[nbnivo]],scannf=FALSE)
    nb.v.p <- length(interres$eig)
    res.afmh <- dudi.pca(X,col.w=poids[[nbnivo]],scannf=FALSE,nf=nb.v.p,scale=FALSE)

plot(res.afmh$li[, coord], asp=1,cex = 0.8,pch=20,xlim=1.1*c(min(res.afmh$li[,coord[1]]),max(res.afmh$li[,coord[1]])),ylim=1.1*c(min(res.afmh$li[,coord[2]]),max(res.afmh$li[,coord[2]])),xlab=paste("Comp",coord[1]," (",signif(res.afmh$eig[coord[1]]*100/sum(res.afmh$eig),4),"%)",sep=""),ylab=paste("Comp",coord[2]," (",signif(res.afmh$eig[coord[2]]*100/sum(res.afmh$eig),4),"%)",sep=""))
abline(v=0,lty=2)
abline(h=0,lty=2)
text(res.afmh$li[,coord[1]],res.afmh$li[,coord[2]],rownames(res.afmh$li), pos = 4,offset = 0.2)
title(main = "Individual map")
get(getOption("device"))()
#senso.corcircle(res.afmh$co[, coord], fullcircle = TRUE)
senso.corcircle(interres$co[, coord], fullcircle = TRUE)
title(paste("Col projection: comp",coord[1]," (",signif(res.afmh$eig[coord[1]]*100/sum(res.afmh$eig),4),"%) - ", "comp", coord[2]," (",signif(res.afmh$eig[coord[2]]*100/sum(res.afmh$eig),4),"%)",sep=""))

    part1 <- list()
    for (h in 1:nbnivo) {
    nbgroup <- length(H[[h]])
    part2 <- list()
    for (g in 1:nbgroup) {
        formule <- matrix(0,dim(X)[1],nb.v.p)
        formule <- (as.matrix(res1[[h]][[g]])) %*% diag(poids[[nbnivo]]) %*% t(X) %*% diag(rep((1/nbind),nbind))
        formule <- formule %*% as.matrix(res.afmh$li) %*% diag(1/res.afmh$eig[1:nb.v.p]) * dilat[[h]][g]
        name1 <- row.names(X)
        name2 <- rep(h,nbind)
        name3 <- rep(g,nbind)
        name <- paste(name1,name2,name3,sep=".")
        namecol1 <- seq(1,nb.v.p)
        namecol <- paste("F",namecol1,sep="")
        formule <- matrix(formule,dim(formule)[1],dim(formule)[2],dimnames=list(name,namecol))
        part2[[g]] <- formule
    }
    part1[[h]] <- part2
    }
    part1[[nbnivo+1]]=res.afmh
    return(part1)
}
###############################################################"

###############################################################"
partial.tab.pour.plot <- function(P,H,coord=c(1,2)) {
    ptpp1 <- list()
    nbnivo <- length(H)
    for (h in 2:nbnivo) {
    x <- 1
    ptpp2 <- list()
    nbgroup <- length(H[[h]])
    for (g in 1:nbgroup) {
        nb <- H[[h]][g]
            ptpp2[[g]]<-P[[h-1]][[x]][,coord]
        if (nb > 1) {
        for (n in 2:nb) {
        ptpp2[[g]] <- rbind(ptpp2[[g]],P[[h-1]][[n+x-1]][,coord])
        }
        }
        x <- x + nb
    }
    ptpp1[[h-1]] <- ptpp2
    }
    ptpp2 <- list()
    nb <- length(H[[nbnivo]])
    ptpp2[[1]] <- P[[nbnivo]][[1]][,coord]
    for (n in 2:nb) {
    ptpp2[[1]] <- rbind(ptpp2[[1]],P[[nbnivo]][[n]][,coord])
    }
    ptpp1[[nbnivo]] <- ptpp2
    return(ptpp1)
}
###############################################################"

###############################################################"
"plot.partial.ind" <- function(coord, H, numb = 6,lab) {

    nbpart1 <- length(coord[[1]])
    inter <- coord[[1]][[1]]
    for(i in 2:nbpart1) {
      inter <- rbind(inter,coord[[1]][[i]])
    }
    xy <- maxmin(inter,H)
    nbnivo <- length(coord)
    nbgroup <- H[[2]][1]
    nbind <- dim(coord[[1]][[1]])[1]/nbgroup
    mult <- nbind / numb
    if ((nbind / numb)!=as.integer(nbind / numb))   mult <- as.integer(mult) + 1

    for (m in 1:mult) {
          par(mfrow = c(2,(numb/2)))
          for (ind in 1:numb) {
        indd <- (m-1)*numb+ind
        if (indd <= nbind) {
            nbgroup <- H[[2]][1]
            nbind <- dim(coord[[1]][[1]])[1]/nbgroup
            fac1 <- factor(rep(1:nbind,nbgroup))
                col_or <- rep(0,length(levels(fac1)))
                col_or[indd] <- 1
            if (((xy[indd,1]>0) & (xy[indd,2]<0)) | ((xy[indd,1]<0) & (xy[indd,2]>0))) {
            ech <- (abs(xy[indd,1])+abs(xy[indd,2]))/10
            }
            else {
                ech <- (xy[indd,2] - xy[indd,1])/10
            }
            if (((xy[indd,3]>0) & (xy[indd,4]<0)) | ((xy[indd,3]<0) & (xy[indd,4]>0))) {
            ech <- (abs(xy[indd,3])+abs(xy[indd,4]))/10
            }
            else {
                ech <- (xy[indd,4] - xy[indd,3])/10
            }
            ech <- abs(ech) 
            x <- c((xy[indd,1]-ech),(xy[indd,2]+ech))
            y <- c((xy[indd,3]-ech),(xy[indd,4]+ech))
                senso.class(coord[[1]][[1]],fac1, xlim=x, ylim=y, cpoint = 0, clabel = 0, col = col_or)
            for (k in 1:(nbnivo-1)) {
                for (j in 1:length(H[[k+1]])) {
                    nbgroup <- H[[k+1]][j]
                    if(nbgroup>1) {
                nbind <- dim(coord[[k]][[j]])[1]/nbgroup
                fac1 <- factor(rep(1:nbind,nbgroup))
                col_or <- rep(0,length(levels(fac1)))
                col_or[indd] <- k
                senso.class(coord[[k]][[j]], fac1, cpoint = 0, clabel = 0,col=col_or,add.plot=TRUE)
                    }
                }
                }
            nbgroup <- length(H[[nbnivo]])
            nbind <- dim(coord[[nbnivo]][[1]])[1]/nbgroup
            fac1 <- factor(rep(1:nbind,nbgroup))
            col_or <- rep(0,length(levels(fac1)))
            col_or[indd] <- nbnivo
            senso.class(coord[[nbnivo]][[1]], fac1, cpoint = 0, clabel = 0, col = col_or, add.plot=TRUE, sub = lab[indd],csub = 1.5)
            }
      }
      if (m < mult)     windows()
    }
}
###############################################################"

###############################################################"
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
###############################################################"

###############################################################"
htabdes <- function(X) {
    nbnivh <- length(X)
    nbvarh <- X
    if (nbnivh>1) {
      for (i in 2:nbnivh) {
        for (j in 1:length(X[[i]])) {    
        nbvarh[[i]][j] <- 0
        if (j==1) {
          for (k in 1:X[[i]][1]) nbvarh[[i]][j] <- nbvarh[[i]][j]+nbvarh[[i-1]][k]
        }
        else {         
            a <- 0     
            b <- 0      
            for (n in 1:(j-1)) a <- a+X[[i]][n]
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
###############################################################"

###############################################################"
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
###############################################################"

"plot.partial" <-
function(coord.partial, H , clab.part = 0, clab.first.part = 0 ,coord=c(1,2),lab){
    nbpart1 <- length(coord.partial[[1]])
    inter <- coord.partial[[1]][[1]]
    for(i in 2:nbpart1) {
      inter <- rbind(inter,coord.partial[[1]][[i]])
    }
        ##
    xmin <- min(inter[,1])-0.2
    xmax <- max(inter[,1])+0.2
    ymin <- min(inter[,2])-0.2
    ymax <- max(inter[,2])+0.2
    x <- c(xmin,xmax)
    y <- c(ymin,ymax)
    nbnivo <- length(coord.partial)
    nbgroup <- H[[2]][1]
    nbind <- dim(coord.partial[[1]][[1]])[1]/nbgroup
    fac1 <- factor(rep(1:nbind,nbgroup))
    senso.class(coord.partial[[1]][[1]],fac1,xlim=x,ylim=y, cpoint = clab.part , clabel = 0,col=rep(1,length(levels(fac1))))
    for (k in 1:(nbnivo-1)) {
      for (j in 1:length(H[[k+1]])) {
            nbgroup <- H[[k+1]][j]
            if(nbgroup>1) {
            nbind <- dim(coord.partial[[k]][[j]])[1]/nbgroup
            fac1 <- factor(rep(1:nbind,nbgroup))
            senso.class(coord.partial[[k]][[j]],fac1, cpoint = clab.part , clabel = 0,col=rep(k,length(levels(fac1))),add.plot=TRUE)
              }
          }
    }
    nbgroup <- length(H[[nbnivo]])
    nbind <- dim(coord.partial[[nbnivo]][[1]])[1]/nbgroup
    fac1 <- factor(rep(lab,nbgroup))
    senso.class(coord.partial[[nbnivo]][[1]],fac1, cpoint = clab.first.part , clabel=0.8, col.moy = rep(nbnivo,length(levels(fac1))),col=rep(nbnivo,length(levels(fac1))),add.plot=TRUE)
    title(sub = "Superimposed representation of the partial clouds" , cex.sub = 1.2, col.sub = "steelblue4", font.sub=3, adj = 0.5, line = 3)
    title(sub = paste("(","comp",coord[1],"-","comp",coord[2],")"), cex.sub = 1.2, col.sub = "steelblue4", font.sub=3, adj = 0.5, line = 4)
}

###############################################################"


###############################################################"
# Main program
###############################################################"

    coord.partial <- partial.ind.tab(X = X, H = H, coord = coord)
    if (length(H)>1){
      tab.coord.partial <- partial.tab.pour.plot(coord.partial, H = H, coord = coord)
      get(getOption("device"))()
      plot.partial(tab.coord.partial, H = H, lab=row.names(X),coord=coord)
      get(getOption("device"))()
      plot.partial.ind(tab.coord.partial ,H = H, num = num,lab=row.names(X))
    }
    res=list()
    res$hmfa=coord.partial[[length(coord.partial)]]
    res$coord.partial=coord.partial[1:(length(coord.partial)-1)]
#    res$tab.coord=tab.coord.partial
    return(res)
}
