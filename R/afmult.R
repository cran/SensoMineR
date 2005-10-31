"afmult" <-function (X, option = c("lambda1", "inertia", "uniform", "internal"),coord=c(1,2),scannf = TRUE, nf = 3,cex=0.8,col = "steelblue4", font = 2, clabel = 0.8,scale.unit=TRUE) {
    if (!inherits(X, "ktab"))   stop("object 'ktab' expected")
    
    
mat = X[[1]]
for (j in 2:length(X$blo)) mat=cbind(mat,X[[j]])
if (any(is.na(mat))){
  for (j in 1:ncol(mat)) mat[,j] <- replace(mat[,j],is.na(mat[,j]),mean(mat[,j],na.rm=TRUE))
}
mat=as.matrix(mat)
un<-rep(1,nrow(mat))
g<-t(mat)%*%diag(X$lw)%*%un
XX<-mat-un%*%t(g)
if (scale.unit==TRUE) {
  InvV<-diag(1/sqrt(diag(t(XX)%*%diag(X$lw)%*%XX)))
  XX<-XX%*%InvV   
}
XX=as.data.frame(XX)
X <- ktab.data.frame(XX,blocks=X$blo,colnames=col.names.ktab(X),rownames=row.names.ktab(X),tabnames=tab.names(X),w.row=X$lw)
    
    if (option[1] == "internal") {
        if (is.null(X$tabw)) {
            warning("Internal weights not found: uniform weigths are used")
            option <- "uniform"
        }
    }
    lw <- X$lw
    cw <- X$cw
    Y <- X
    sepan <- sepan(X, nf = 4)
    nbloc <- length(sepan$blo)
    indicablo <- factor(rep(1:nbloc, sepan$blo))
    rank.fac <- factor(rep(1:nbloc, sepan$rank))
    ncw <- NULL
    tab.names <- names(X)[1:nbloc]
    auxinames <- ktab.util.names(X)
    if (option[1] == "lambda1") {
        for (i in 1:nbloc) {
            ncw <- c(ncw, rep(1/sepan$Eig[rank.fac == i][1], sepan$blo[i]))
        }
    }
    else if (option[1] == "inertia") {
        for (i in 1:nbloc) {
            ncw <- c(ncw, rep(1/sum(sepan$Eig[rank.fac == i]), sepan$blo[i]))
        }
    }
    else if (option[1] == "uniform")    ncw <- rep(1, sum(sepan$blo))
    else if (option[1] == "internal")    ncw <- rep(X$tabw, sepan$blo)
    else stop("unknown option")
    ncw <- cw * ncw
    tab <- X[[1]]
    for (i in 2:nbloc) {
        tab <- cbind.data.frame(tab, X[[i]])
    }
    names(tab) <- auxinames$col
    anaco <- as.dudi(tab, col.w = ncw, row.w = lw, nf = nf, scannf = scannf, call = match.call(), type = "mfa")
    nf <- anaco$nf
    afm <- list()
    afm$tab.names <- names(X)[1:nbloc]
    afm$blo <- X$blo
    afm$TL <- X$TL
    afm$TC <- X$TC
    afm$T4 <- X$T4
    afm$tab <- anaco$tab
    afm$eig <- anaco$eig
    afm$rank <- anaco$rank
    afm$li <- anaco$li
    afm$l1 <- anaco$l1
    afm$nf <- anaco$nf
    afm$lw <- anaco$lw
    afm$cw <- anaco$cw
    afm$co <- anaco$co
    afm$c1 <- anaco$c1
#    afm$seb <- anaco
    projiner <- function(xk, qk, d, z) {
        w7 <- t(as.matrix(xk) * d) %*% as.matrix(z)
        iner <- apply(w7 * w7 * qk, 2, sum)
        return(iner)
    }
    link <- matrix(0, nbloc, nf)
    for (k in 1:nbloc) {
        xk <- X[[k]]
        q <- ncw[indicablo == k]
        link[k, ] <- projiner(xk, q, lw, anaco$l1)
    }
    link <- as.data.frame(link)
    names(link) <- paste("Comp", 1:nf, sep = "")
    row.names(link) <- tab.names
    afm$link <- link
    w <- matrix(0, nbloc * 4, nf)
    i1 <- 0
    i2 <- 0
    matl1 <- as.matrix(afm$l1)
    matli <- as.matrix(afm$li)
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + 4
        tab <- as.matrix(sepan$L1[sepan$TL[, 1] == k, ])
        if (ncol(tab) > 4) 
            tab <- tab[, 1:4]
        if (ncol(tab) < 4)   tab <- cbind(tab, matrix(0, nrow(tab), 4 - ncol(tab)))
        tab <- t(tab * lw) %*% matl1
        for (i in 1:min(nf, 4)) {
            if (tab[i, i] < 0) {
                for (j in 1:nf) tab[i, j] <- -tab[i, j]
            }
        }
        w[i1:i2, ] <- tab
    }
    w <- data.frame(w)
    names(w) <- paste("Comp", 1:nf, sep = "")
    row.names(w) <- auxinames$tab
    afm$T4comp <- w
    poids <- diag(1/afm$eig[1:nf])
    w <- matrix(0, nrow(sepan$TL), ncol = nf)
    i1 <- 0
    i2 <- 0
    for (k in 1:nbloc) {
        i1 <- i2 + 1
        i2 <- i2 + length(lw)
        qk <- ncw[indicablo == k]
        xk <- as.matrix(X[[k]])
        w[i1:i2, ] <- (xk %*% (qk * t(xk))) %*% (diag(afm$lw) %*% matli %*% poids * nbloc)
      }
    w <- data.frame(w)
    row.names(w) <- auxinames$row
    names(w) <- paste("Fac", 1:nf, sep = "")
    afm$lisup <- w
    afm$tabw <- X$tabw
    afm$call <- match.call()
#########################################################################
#scatterutil.eigen(anaco$eig, sub = "", wsel = coord)
#title(main = "Eigenvalues", cex.main = cex, font.main = font,col.main = col, adj = 1)
#get(getOption("device"))()
#par()

plot(afm$li[, coord], asp=1,cex = 0.8,pch=20,xlab=paste("Comp",coord[1]," (",signif(anaco$eig[coord[1]]*100/sum(anaco$eig),4),"%)",sep=""),ylab=paste("Comp",coord[2]," (",signif(anaco$eig[coord[2]]*100/sum(anaco$eig),4),"%)",sep=""))
abline(v=0,lty=2)
abline(h=0,lty=2)
text(afm$li[,coord[1]],afm$li[,coord[2]],rownames(afm$li), pos = 4,offset = 0.2,col=col)
title(sub = "Row projection", cex.sub = cex, font.sub = font,col.sub = col, adj = 0, line = 3.8)
get(getOption("device"))()
coolig <- w[, coord]
senso.class(coolig, fac = as.factor(afm$TL[, 2]),clabel = 0,box=FALSE)
text(afm$li[,coord[1]],afm$li[,coord[2]],rownames(afm$li), pos = 4,offset = 0.2,col=col)
title(sub = paste("Row projection (and their partial representations): comp",coord[1]," (",signif(anaco$eig[coord[1]]*100/sum(anaco$eig),4),"%) - ", "comp", coord[2]," (",signif(anaco$eig[coord[2]]*100/sum(anaco$eig),4),"%)",sep=""), cex.sub = cex, font.sub = font,col.sub = col, adj = 0.5, line = 3.8)
get(getOption("device"))()
senso.corcircle(afm$co[, coord], fullcircle = TRUE,clabel = clabel)
title(sub = paste("Col projection: comp",coord[1]," (",signif(anaco$eig[coord[1]]*100/sum(anaco$eig),4),"%) - ", "comp", coord[2]," (",signif(anaco$eig[coord[2]]*100/sum(anaco$eig),4),"%)",sep=""), cex.sub = cex, font.sub = font,col.sub = col, adj = 0.5, line = 3.8)
get(getOption("device"))()
par(mar = c(5, 4, 2, 2))
plot(afm$link[, coord], xlim = c(0, 1), ylim = c(0,1), asp = 1, pch = 3,xlab=paste("Comp",coord[1]," (",signif(anaco$eig[coord[1]]*100/sum(anaco$eig),4),"%)",sep=""),ylab=paste("Comp",coord[2]," (",signif(anaco$eig[coord[2]]*100/sum(anaco$eig),4),"%)",sep=""))
scatterutil.grid(0)
title(sub = "Link", cex.sub = cex, font.sub = font,col.sub = col, adj = 0, line = 3.6)
senso.scatterutil.eti(afm$link[, coord[1]], afm$link[,coord[2]], label = row.names(afm$link), clabel = clabel,pos = 3, font = 1)
get(getOption("device"))()
#par(mar = c(5, 4, 4, 2))
senso.corcircle(afm$T4comp[afm$T4[, 2] == 1,coord], fullcircle = TRUE, box = FALSE, clabel = clabel)
for (j in 2:max(coord)) senso.corcircle(afm$T4comp[afm$T4[, 2] == j,coord], fullcircle = TRUE, box = FALSE, clabel = clabel,add.plot=TRUE)
title(sub = paste("Component projection: comp",coord[1]," (",signif(anaco$eig[coord[1]]*100/sum(anaco$eig),4),"%) - ", "comp", coord[2]," (",signif(anaco$eig[coord[2]]*100/sum(anaco$eig),4),"%)",sep=""), cex.sub = cex,font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
#########################################################################
    class(afm) <- c("mfa", "list")
    return(afm)
    return(rank.fac)
}
