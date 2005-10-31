carto <- function (Mat, MatH, option = c("pca", "manual","afmult", "statis","gpa"), level = 0, regmod = 1, graphic = FALSE,
    coord = c(1, 2), asp = 1, cex = 1.3, col = "steelblue4", font = 2, clabel = 0.8,label.j=FALSE,resolution=200,nb.clusters=0)
{
get(getOption("device"))()
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    predire <- function(n1, n2, coeff) {
        coeff[1] + coeff[2] * n1 + coeff[3] * n2 + coeff[4] *
            n1 * n1 + coeff[5] * n2 * n2 + coeff[6] * n1 * n2
    }
    if (!is.data.frame(MatH))
        stop("Non convenient selection for MatH")
    if (option[1] == "pca") {
        if (!is.data.frame(Mat))
            stop("Non convenient selection for Mat")
        interpca <- dudi.pca(Mat, scannf = FALSE)
        nb.v.p <- length(interpca$eig)
        pcamap <- dudi.pca(Mat, scannf = FALSE, nf = nb.v.p)
        if (!graphic) {
            par(mfrow = c(2, 1))
            senso.label(pcamap$li[, coord], clabel = 0, cpoint = 0.8)
            add.scatter.eig(pcamap$eig, pcamap$nf, coord[1],
                coord[2], posi = "top", ratio = 1/5)
            text(pcamap$li[, coord[1]], pcamap$li[, coord[2]],
                labels = row.names(Mat), cex = clabel, pos = 4,
                offset = 0.2)
            title(sub = paste("Individuals factor map (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.8)
            senso.corcircle(pcamap$co[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Correlation circle (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.8)
            get(getOption("device"))()
        }
        else if (graphic) {
            scatterutil.eigen(pcamap$eig, sub = "", wsel = coord)
            title(main = "Eigenvalues", cex.main = cex, font.main = font,
                col.main = col, adj = 1)
            get(getOption("device"))()
            par(mar = c(0, 0, 2, 0))
            senso.label(pcamap$li[, coord], clabel = 0, cpoint = 0.8)
            text(pcamap$li[, coord[1]], pcamap$li[, coord[2]],
                labels = row.names(Mat), cex = clabel, pos = 4,
                offset = 0.2)
            title(main = paste("Individuals factor map (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.main = cex,
                font.main = font, col.main = col, adj = 0.5)
            get(getOption("device"))()
            senso.corcircle(pcamap$co[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Correlation circle (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
            get(getOption("device"))()
        }
        matrice <- cbind(row.names(MatH), pcamap$li[, coord],MatH)
  aa=cor(matrice[,4:ncol(matrice)],matrice[,2:3])
  bb=cor(Mat,matrice[,2:3])
  hc <- hclust(dist(t(MatH)),method="ward")
  plot(as.dendrogram(hc),main="Cluster Dendrogram",xlab="Panelists",leaflab="none") 
  get(getOption("device"))()
  if (nb.clusters==0){
    classif=hopach(t(MatH),d="cor",K=10,mss="mean")
    nb.clusters=classif$clustering$k
  }
  aux=pam(t(MatH),k=nb.clusters)$clustering
  mat <- matrix(0,nb.clusters,nrow(MatH))
  dimnames(mat) <- list(1:nb.clusters,rownames(MatH))
  for (i in 1:nb.clusters){
    mat[i,] <- apply(t(MatH[,aux==i]),2,mean)
    rownames(mat)[i] <- paste("cluster",i)
  }  
  ab=cor(t(mat),matrice[,2:3])
  senso.corcircle(bb, fullcircle=TRUE,col="red")
  if (label.j==FALSE) senso.corcircle(aa, fullcircle=TRUE,label=NULL,lty=2,add=TRUE)
  if (label.j==TRUE) senso.corcircle(aa, fullcircle=TRUE,lty=2,add=TRUE)
  senso.corcircle(ab, fullcircle=TRUE,col="blue",add=TRUE)
  title(main = paste("Correlation circle (comp",coord[1],"-","comp",coord[2],")"))
  get(getOption("device"))()
    }
    else if (option[1] == "manual") {
        matrice <- cbind(row.names(MatH), Mat,MatH)
   }
    else if (option[1] == "afmult") {
        if (!inherits(Mat, "ktab"))
            stop("Object of type 'ktab' expected")
        intermfa <- afmult(Mat, scannf = FALSE)
        nb.v.p <- length(intermfa$eig)
        mfamap <- afmult(Mat, scannf = FALSE, nf = nb.v.p)
        if (!graphic) {
            par(mfrow = c(2, 2))
            coolig <- mfamap$lisup[, coord]
            senso.class(coolig, fac = as.factor(mfamap$TL[, 2]),clabel = 0)
            text(mfamap$li[,coord[1]],mfamap$li[,coord[2]],rownames(mfamap$li), pos = 4,offset = 0.2,col=col)
            add.scatter.eig(mfamap$eig, mfamap$nf, coord[1],
                coord[2], posi = "top", ratio = 1/5)
            title(sub = paste("Row projection (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.8)
            senso.corcircle(mfamap$co[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Col projection (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.8)
            par(mar = c(5, 4, 2, 2))
            plot(mfamap$link[, coord], type = "s", xlim = c(0,
                1), ylim = c(0, 1), asp = 1, pch = 3)
            scatterutil.grid(0)
            title(sub = "Link", cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0, line = 3.6)
            senso.scatterutil.eti(mfamap$link[, coord[1]], mfamap$link[,
                coord[2]], label = row.names(mfamap$link), clabel = clabel,
                pos = 3, font = 1)
            par(mar = c(5, 4, 4, 2))
            senso.corcircle(mfamap$T4comp[mfamap$T4[, 2] == 1,
                coord], fullcircle = TRUE, box = TRUE, clabel = clabel)
            title(sub = paste("Component projection (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
        }
        else if (graphic) {
            scatterutil.eigen(mfamap$eig, sub = "", wsel = coord)
            title(main = "Eigenvalues", cex.main = cex, font.main = font,
                col.main = col, adj = 1)
            get(getOption("device"))()
            par()
            coolig <- mfamap$lisup[, coord]
            senso.class(coolig, fac = as.factor(mfamap$TL[, 2]),clabel = 0)
            text(mfamap$li[,coord[1]],mfamap$li[,coord[2]],rownames(mfamap$li), pos = 4,offset = 0.2,col=col)
            title(sub = paste("Row projection (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.8)
            get(getOption("device"))()
            senso.corcircle(mfamap$co[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Col projection (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.8)
            get(getOption("device"))()
            par(mar = c(5, 4, 2, 2))
            plot(mfamap$link[, coord], xlim = c(0, 1), ylim = c(0,
                1), asp = 1, pch = 3)
            scatterutil.grid(0)
            title(sub = "Link", cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0, line = 3.6)
            senso.scatterutil.eti(mfamap$link[, coord[1]], mfamap$link[,
                coord[2]], label = row.names(mfamap$link), clabel = clabel,
                pos = 3, font = 1)
            get(getOption("device"))()
            par(mar = c(5, 4, 4, 2))
            senso.corcircle(mfamap$T4comp[mfamap$T4[, 2] == 1,
                coord], fullcircle = TRUE, box = TRUE, clabel = clabel)
            title(sub = paste("Component projection (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
        }
        get(getOption("device"))()
        matrice <- cbind(row.names(MatH), mfamap$li[, coord],
            MatH)
    }
    else if (option[1] == "statis") {
        if (!inherits(Mat, "ktab"))
            stop("Object of type 'ktab' expected")
        interstatis <- statis(Mat, scannf = FALSE)
        nb.v.p <- length(interstatis$eig)
        statismap <- statis(Mat, scannf = FALSE, nf = nb.v.p)
        if (!graphic) {
            par(mfrow = c(2, 2))
            par(mar = c(5, 4, 4, 2))
            senso.corcircle(statismap$RV.coo[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Interstructure (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.6)
            par(mar = c(5.5, 0, 0, 0))
            senso.label(statismap$C.li[, coord], clabel = 0, cpoint = 0.8)
            text(statismap$C.li[, coord[1]], statismap$C.li[,
                coord[2]], labels = row.names(Mat), cex = clabel,
                pos = 1, offset = 0.2)
            add.scatter.eig(statismap$RV.eig, length(statismap$RV.eig),
                coord[1], coord[2], posi = "top", ratio = 1/4)
            title(sub = paste("Compromise (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5)
            par(mar = c(5, 4, 4, 2))
            senso.corcircle(statismap$C.T4[statismap$T4[, 2] ==
                1, ], fullcircle = TRUE, clabel = clabel, box = TRUE)
            title(sub = paste("Component projection (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
            par(mar = c(5.5, 4, 1.5, 1.5))
            plot(statismap$RV.tabw, statismap$cos2, xlab = "Tables weights",
                ylab = "Cos 2", xlim = c(0, 1), ylim = c(0, 1),
                asp = 1, pch = 3, frame.plot = TRUE)
            scatterutil.grid(0)
            title(sub = "Typological value", cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0)
            senso.scatterutil.eti(statismap$RV.tabw, statismap$cos2,
                label = statismap$tab.names, clabel = clabel,
                pos = 3, font = 1)
            get(getOption("device"))()
        }
        else if (graphic) {
            scatterutil.eigen(statismap$RV.eig, sub = "", wsel = coord)
            title(main = "Eigenvalues", cex.main = cex, font.main = font,
                col.main = col, adj = 1)
            get(getOption("device"))()
            senso.corcircle(statismap$RV.coo[, coord], fullcircle = TRUE,
                clabel = clabel)
            title(sub = paste("Interstructure (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5, line = 3.6)
            get(getOption("device"))()
            par(mar = c(5.5, 0, 0, 0))
            senso.label(statismap$C.li[, coord], clabel = 0, cpoint = 0.8)
            text(statismap$C.li[, coord[1]], statismap$C.li[,
                coord[2]], labels = row.names(Mat), cex = clabel,
                pos = 1, offset = 0.2)
            title(sub = paste("Compromise (", "comp", coord[1],
                "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0.5)
            get(getOption("device"))()
            par(mar = c(5, 4, 4, 2))
            senso.corcircle(statismap$C.T4[statismap$T4[, 2] ==
                1, ], fullcircle = TRUE, clabel = clabel, box = TRUE)
            title(sub = paste("Component projection (", "comp",
                coord[1], "-", "comp", coord[2], ")"), cex.sub = cex,
                font.sub = font, col.sub = col, adj = 0.5, line = 3.6)
            get(getOption("device"))()
            par(mar = c(5.5, 4, 1.5, 1.5))
            plot(statismap$RV.tabw, statismap$cos2, xlab = "Tables weights",
                ylab = "Cos 2", xlim = c(0, 1), ylim = c(0, 1),
                asp = 1, pch = 3, frame.plot = TRUE)
            scatterutil.grid(0)
            title(sub = "Typological value", cex.sub = cex, font.sub = font,
                col.sub = col, adj = 0)
            senso.scatterutil.eti(statismap$RV.tabw, statismap$cos2,
                label = statismap$tab.names, clabel = clabel,
                pos = 3, font = 1)
            get(getOption("device"))()
        }
        matrice <- cbind(row.names(MatH), statismap$C.li[, coord],
            MatH)
    }
    else if (option[1] == "gpa") {
        if (!inherits(Mat, "ktab"))
            stop("Object of type 'ktab' expected")
        resultats <- gpa(Mat)
        par()
        res <- resultats$Xfin[, , 1]
        for (i in 2:dim(resultats$Xfin)[3]) {
            res <- rbind(res, resultats$Xfin[, , i])
        }
        coo.part <- data.frame(res[, coord])
        row.names(coo.part) <- paste(rep(1:dim(resultats$Xfin)[1],
            dim(resultats$Xfin)[3]), rep(1:dim(resultats$Xfin)[3],
            each = dim(resultats$Xfin)[1]), sep = ".")
        fact.gpa <- factor(rep((1:dim(resultats$Xfin)[1]), dim(resultats$Xfin)[3]))
        resultats$consensus <- data.frame(resultats$consensus)
        lab.gpa <- row.names(resultats$consensus)
        senso.class(coo.part, fac = fact.gpa, clabel = clabel)
        if (!graphic) {
            add.scatter.eig(diag(resultats$gama), length(diag(resultats$gama)),
                coord[1], coord[2], posi = "top", ratio = 1/5)
        }
        title(sub = paste("Row projection (", "comp", coord[1],
            "-", "comp", coord[2], ")"), cex.sub = cex, font.sub = font,
            col.sub = col, adj = 0.5, line = 3.8)
        get(getOption("device"))()
        matrice <- cbind(row.names(MatH), resultats$consensus[,
            coord], MatH)
    }

  if (option[1]!="dudi.pca"){
    hc <- hclust(dist(t(MatH)),method="ward")
    plot(as.dendrogram(hc),main="Cluster Dendrogram",xlab="Panelists",leaflab="none") 
    get(getOption("device"))()
    if (nb.clusters==0){
      classif=hopach(t(MatH),d="cor",K=10,mss="mean")
      nb.clusters=classif$clustering$k
    }
    aux=pam(t(MatH),k=nb.clusters)$clustering
    mat <- matrix(0,nb.clusters,nrow(MatH))
    dimnames(mat) <- list(1:nb.clusters,rownames(MatH))
    for (i in 1:nb.clusters){
      mat[i,] <- apply(t(MatH[,aux==i]),2,mean)
      rownames(mat)[i] <- paste("cluster",i)
    }  
    ab=cor(t(mat),matrice[,2:3])
    aa=cor(matrice[,4:ncol(matrice)],matrice[,2:3])
    if (label.j==FALSE) senso.corcircle(aa, fullcircle=TRUE,lty=2,label=NULL)
    if (label.j==TRUE) senso.corcircle(aa, fullcircle=TRUE,lty=2)
     senso.corcircle(ab, fullcircle=TRUE,col="blue",add=TRUE)
    title(main = paste("Correlation circle (comp",coord[1],"-","comp",coord[2],")"))
    get(getOption("device"))()
  }
  
    matrice[, 4:ncol(matrice)] <- scale(matrice[, 4:ncol(matrice)],
        center =TRUE, scale =FALSE)[, ]
    nbconso <- ncol(matrice) - 3
    x1 <- matrice[, 2]
    x2 <- matrice[, 3]
    x12 <- scale(x1, center =TRUE, scale =FALSE)[, ]^2
    x22 <- scale(x2, center =TRUE, scale =FALSE)[, ]^2
    x12plusx22 <- x12 + x22
    x3 <- scale(x1, center =TRUE, scale =FALSE)[, ] * scale(x2, center =TRUE,
        scale =FALSE)[, ]
    XX <- cbind(x1, x2, x12, x22, x3)
    etendue.x1 <- diff(range(x1))
    etendue.x2 <- diff(range(x2))
    pas <- max(etendue.x1, etendue.x2)/resolution
    f1 <- seq((min(x1)-etendue.x1*0.05),(max(x1)+etendue.x1*0.05), pas)
    f2 <- seq((min(x2)-etendue.x2*0.05),(max(x2)+etendue.x2*0.05), pas)
    depasse <- matrix(0, nrow = length(f1), ncol = length(f2))
    abscis <- NULL
    ordon <- NULL
    for (i in 1:nbconso) {
        if (regmod == 1)
            coeff <- lm(matrice[, i + 3] ~ XX[, 1] + XX[, 2] +
                XX[, 3] + XX[, 4] + XX[, 5])$coef
        if (regmod == 2) {
            coeff <- lm(matrice[, i + 3] ~ XX[, 1] + XX[, 2])$coef
            coeff <- c(coeff, 0, 0, 0)
        }
        if (regmod == 3) {
            coeff <- lm(matrice[, i + 3] ~ x1 + x2 + x12plusx22)$coef
            coeff <- c(coeff, coeff[4], 0)
        }
        if (regmod == 4) {
            coeff <- lm(matrice[, i + 3] ~ XX[, 1] + XX[, 2] +
                XX[, 3] + XX[, 4])$coef
            coeff <- c(coeff, 0)
        }
        predites <- outer(f1, f2, predire, coeff)
        predites <- (predites - mean(predites, na.rm =TRUE))/sqrt(var(as.vector(predites),
            na.rm =TRUE))
        depasse <- depasse + matrix(as.numeric(predites > level),
            nrow = length(f1), ncol = length(f2))
        abscis <- c(abscis, f1[rev(order(predites))[1] - length(f1) *
            as.integer((rev(order(predites))[1] - 0.5)/length(f1))])
        ordon <- c(ordon, f2[as.integer(1 + (rev(order(predites))[1] -
            0.5)/length(f1))])
    }
    depasse <- round(depasse/nbconso * 100)
    dimnames(depasse) <- list(as.character(f1), as.character(f2))
    image(f1, f2, depasse, col = terrain.colors(200), xlab = "Comp 1",
        ylab = "Comp 2", main = "Preference mapping", font.main = font,
        col.main = col, cex.main = cex, asp = asp)
    contour(f1, f2, depasse, nlevels = 9, levels = c(20, 30,
        40, 50, 60, 70, 80, 90, 95), add =TRUE, labex = 0)
    for (i in 1:nrow(matrice)) {
        points(matrice[i, 2], matrice[i, 3])
        text(matrice[i, 2] + 0.01, matrice[i, 3] + 0.01, matrice[i,
            1])
    }
    points(abscis, ordon, pch = 20)
}
