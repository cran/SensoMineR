optimaldesign = function (nbPanelist, nbProd, nbProdByPanelist, nbPanelistMin = nbPanelist,
    ordre = TRUE, weight = 0.5, graine = Sys.time(), nbDesignProd = 10,
    nbDesignOrdre = 50, matEssImp = NA){

    condenseProd <- function(nbPanelist, nbProd, nbEssai = NULL) {
        colJuge <- gl(nbPanelist, nbProd)
        colProd <- as.factor(seq(nbProd))
        if (is.null(nbEssai) | nbEssai > (nbPanelist * nbProd))
            nbEssai <- nbPanelist * nbProd
        res <- cbind.data.frame(colJuge, colProd)
        return(res[1:nbEssai, ])
    }
    matEffet <- function(nbPanelist, nbProd, matCond) {
        res <- model.matrix(~colProd + colJuge, data = matCond,
            contrasts = list(colProd = "contr.sum", colJuge = "contr.sum"))
        partieProd <- 2:nbProd
        partieJuge <- c((nbProd + 1):(nbProd + nbPanelist - 1))
        res <- res[, c(partieJuge, partieProd, 1)]
        return(res)
    }
    genJPAlea <- function(nbEssaiImp, indEssImp, nbEssaiPlan,
        nbEssaiPlanComp) {
        res <- vector(mode = "integer", length = nbEssaiPlanComp)
        dejaTire <- vector(length = nbEssaiPlanComp)
        if (nbEssaiImp > 0) {
            res[1:nbEssaiImp] <- indEssImp[1:nbEssaiImp]
            dejaTire[res[1:nbEssaiImp]] <- TRUE
        }
        nbTirage <- nbEssaiPlan - nbEssaiImp
        i <- 0
        while (i < nbTirage) {
            tirage <- as.integer(1 + runif(1) * nbEssaiPlanComp)
            if (dejaTire[tirage])
                next
            res[nbEssaiImp + 1 + i] <- tirage
            dejaTire[tirage] <- TRUE
            i <- i + 1
        }
        return(res)
    }
    calcEssaisImp <- function(plan) {
        nbLig <- nrow(plan)
        nbCol <- ncol(plan)
        matEssImp <- matrix(nrow = nbLig * nbCol, ncol = 2)
        k <- 0
        for (i in 1:nbLig) {
            for (j in 1:nbCol) {
                if (plan[i, j] == 1) {
                  k <- k + 1
                  matEssImp[k, 1] <- i
                  matEssImp[k, 2] <- j
                }
            }
        }
        return(list(nb = k, valeurs = matEssImp[1:k, ]))
    }
    findEssImp <- function(matEssImp, matCond) {
        nbEssImp = nrow(matEssImp)
        res <- vector(mode = "integer", length = nbEssImp)
        for (ei in 1:nbEssImp) res[ei] <- which(matCond[, 1] ==
            matEssImp[ei, 1] & matCond[, 2] == matEssImp[ei,
            2])
        return(res)
    }
    rang <- function(nbPanelist, nbProd, nbRang, plan, weightRang,
        nbPanelistImp, matJugeImp, nbDesign) {
        genJRAlea <- function(matCond, nbPanelist, nbRang, nbPanelistImp) {
            matAlea <- matCond
            for (i in (nbPanelistImp + 1):nbPanelist) {
                tot <- nbRang
                for (r in 1:(nbRang - 1)) {
                  jal <- as.integer(1 + runif(1) * tot)
                  matAlea[i, r] <- matCond[i, jal]
                  if (jal <= tot - 1) {
                    for (k in jal:(tot - 1)) matCond[i, k] <- matCond[i,
                      k + 1]
                  }
                  tot <- tot - 1
                }
                matAlea[i, nbRang] <- matCond[i, 1]
            }
            return(matAlea)
        }
        condenseRang <- function(nbPanelist, nbProd, nbRang,
            nbPanelistImp, matJugeImp, plan) {
            if (nbPanelistImp != 0) {
                if (nbPanelistImp != nrow(matJugeImp))
                  stop("problem with the judge parameters")
            }
            jugeProCond <- matrix(0, nbPanelist, nbRang)
            if (nbPanelistImp == 0) {
                for (i in 1:nbPanelist) {
                  k <- 1
                  for (j in 1:nbProd) {
                    if (plan[i, j] == 1) {
                      jugeProCond[i, k] <- j
                      k <- k + 1
                    }
                  }
                }
            }
            else {
                jugeProCond[1:nbPanelistImp, ] <- matJugeImp
                for (i in ((nbPanelistImp + 1):nbPanelist)) {
                  k <- 1
                  for (j in 1:nbProd) {
                    if (plan[i, j] == 1) {
                      jugeProCond[i, k] <- j
                      k <- k + 1
                    }
                  }
                }
            }
            return(jugeProCond)
        }
        eam <- function(mat, effTh) {
            ecartUnif <- sum(abs(effTh - mat))
            return(ecartUnif/(nrow(mat) * ncol(mat)))
        }
        majJR <- function(matJR, rang1Max, rang2Max, i) {
            tmp <- matJR[i, rang1Max]
            matJR[i, rang1Max] <- matJR[i, rang2Max]
            matJR[i, rang2Max] <- tmp
            return(matJR)
        }
        majPP <- function(matPP, i, matJR, nbRang, rang1Max,
            rang2Max) {
            if (rang1Max != 1)
                matPP[matJR[i, rang1Max - 1], matJR[i, rang1Max]] <- matPP[matJR[i,
                  rang1Max - 1], matJR[i, rang1Max]] - 1
            if (rang1Max + 1 == rang2Max)
                matPP[matJR[i, rang1Max], matJR[i, rang2Max]] <- matPP[matJR[i,
                  rang1Max], matJR[i, rang2Max]] - 1
            else {
                matPP[matJR[i, rang1Max], matJR[i, rang1Max +
                  1]] <- matPP[matJR[i, rang1Max], matJR[i, rang1Max +
                  1]] - 1
                matPP[matJR[i, rang2Max - 1], matJR[i, rang2Max]] <- matPP[matJR[i,
                  rang2Max - 1], matJR[i, rang2Max]] - 1
            }
            if (rang2Max != nbRang)
                matPP[matJR[i, rang2Max], matJR[i, rang2Max +
                  1]] <- matPP[matJR[i, rang2Max], matJR[i, rang2Max +
                  1]] - 1
            if (rang2Max != nbRang)
                matPP[matJR[i, rang1Max], matJR[i, rang2Max +
                  1]] <- matPP[matJR[i, rang1Max], matJR[i, rang2Max +
                  1]] + 1
            if (rang1Max != 1)
                matPP[matJR[i, rang1Max - 1], matJR[i, rang2Max]] <- matPP[matJR[i,
                  rang1Max - 1], matJR[i, rang2Max]] + 1
            if (rang2Max != rang1Max + 1) {
                matPP[matJR[i, rang2Max], matJR[i, rang1Max +
                  1]] <- matPP[matJR[i, rang2Max], matJR[i, rang1Max +
                  1]] + 1
                matPP[matJR[i, rang2Max - 1], matJR[i, rang1Max]] <- matPP[matJR[i,
                  rang2Max - 1], matJR[i, rang1Max]] + 1
            }
            else matPP[matJR[i, rang2Max], matJR[i, rang1Max]] <- matPP[matJR[i,
                rang2Max], matJR[i, rang1Max]] + 1
            return(matPP)
        }
        majPR <- function(matPR, rang1Max, rang2Max, i, matJR,
            prod1Max, prod2Max) {
            if (matJR[i, rang1Max] < matJR[i, rang2Max]) {
                matPR[prod1Max, rang1Max] <- matPR[prod1Max,
                  rang1Max] - 1
                matPR[prod1Max, rang2Max] <- matPR[prod1Max,
                  rang2Max] + 1
                matPR[prod2Max, rang1Max] <- matPR[prod2Max,
                  rang1Max] + 1
                matPR[prod2Max, rang2Max] <- matPR[prod2Max,
                  rang2Max] - 1
            }
            else {
                matPR[prod1Max, rang1Max] <- matPR[prod1Max,
                  rang1Max] + 1
                matPR[prod1Max, rang2Max] <- matPR[prod1Max,
                  rang2Max] - 1
                matPR[prod2Max, rang1Max] <- matPR[prod2Max,
                  rang1Max] - 1
                matPR[prod2Max, rang2Max] <- matPR[prod2Max,
                  rang2Max] + 1
            }
            return(matPR)
        }
        valOptiPlan <- function(nbPanelist, nbProd, nbRang, nbEssPlan,
            unifPP, unifPR, weightRang, weightSucc, matPR, matPP) {
            optiRang1 <- sum(abs((nbEssPlan/nbRang) - colSums(matPR)))
            optiRang2 <- sum(as.integer(abs(unifPR - matPR)))
            optiRang <- optiRang1 + optiRang2
            totslig <- vector(mode = "integer", length = nbProd)
            for (i in 1:nbProd) {
                for (j in 1:nbProd) {
                  if (i != j)
                    totslig[i] <- totslig[i] + matPP[i, j]
                }
            }
            optiSucc1 <- sum(as.integer(abs(((nbEssPlan - nbPanelist)/nbProd) -
                totslig)))
            optiSucc2 <- (abs(unifPP - matPP))
            diag(optiSucc2) <- 0
            optiSucc2 <- sum(as.integer(optiSucc2))
            optiSucc <- optiSucc1 + optiSucc2
            return(weightRang * optiRang + weightSucc * optiSucc)
        }
        critQualPP <- function(rang1, rang2, nbRang, matPP, unifPP,
            matJR, juge) {
            distuv <- 0
            distup <- 0
            if (rang1 != 1)
                distuv <- distuv + (matPP[matJR[juge, rang1 -
                  1], matJR[juge, rang1]] - unifPP)^2
            if (rang1 + 1 == rang2) {
                distuv <- distuv + (matPP[matJR[juge, rang1],
                  matJR[juge, rang2]] - unifPP)^2
                distuv <- distuv + (unifPP - matPP[matJR[juge,
                  rang2], matJR[juge, rang1]])^2
            }
            else {
                distuv <- distuv + (matPP[matJR[juge, rang1],
                  matJR[juge, rang1 + 1]] - unifPP)^2
                distuv <- distuv + (matPP[matJR[juge, rang2 -
                  1], matJR[juge, rang2]] - unifPP)^2
            }
            if (rang2 != nbRang)
                distuv <- distuv + (matPP[matJR[juge, rang2],
                  matJR[juge, rang2 + 1]] - unifPP)^2
            if (rang1 + 1 != rang2)
                distuv <- distuv + (unifPP - matPP[matJR[juge,
                  rang2 - 1], matJR[juge, rang1]])^2
            if (rang2 != nbRang)
                distuv <- distuv + (unifPP - matPP[matJR[juge,
                  rang1], matJR[juge, rang2 + 1]])^2
            if (rang1 != 1)
                distuv <- distuv + (unifPP - matPP[matJR[juge,
                  rang1 - 1], matJR[juge, rang2]])^2
            if (rang2 != rang1 + 1)
                distuv <- distuv + (unifPP - matPP[matJR[juge,
                  rang2], matJR[juge, rang1 + 1]])^2
            if (rang1 != 1)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang1 - 1], matJR[juge, rang1]] + 1)^2
            if (rang1 + 1 == rang2) {
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang1], matJR[juge, rang2]] + 1)^2
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang2], matJR[juge, rang1]] - 1)^2
            }
            else {
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang1], matJR[juge, rang1 + 1]] + 1)^2
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang2 - 1], matJR[juge, rang2]] + 1)^2
            }
            if (rang2 != nbRang)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang2], matJR[juge, rang2 + 1]] + 1)^2
            if (rang1 + 1 != rang2)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang2 - 1], matJR[juge, rang1]] - 1)^2
            if (rang2 != nbRang)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang1], matJR[juge, rang2 + 1]] - 1)^2
            if (rang1 != 1)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang1 - 1], matJR[juge, rang2]] - 1)^2
            if (rang2 != rang1 + 1)
                distup <- distup + (unifPP - matPP[matJR[juge,
                  rang2], matJR[juge, rang1 + 1]] - 1)^2
            return(distuv - distup)
        }
        critQualPR <- function(prod1, rang1, prod2, rang2, matPR,
            unifPR, matJR, juge) {
            prod1Rang1 <- matPR[prod1, rang1]
            prod1Rang2 <- matPR[prod1, rang2]
            prod2Rang1 <- matPR[prod2, rang1]
            prod2Rang2 <- matPR[prod2, rang2]
            ecav <- (prod1Rang1 - unifPR)^2 + (prod1Rang2 - unifPR)^2 +
                (prod2Rang1 - unifPR)^2 + (prod2Rang2 - unifPR)^2
            if (matJR[juge, rang1] < matJR[juge, rang2]) {
                p1 <- prod1Rang1 - 1
                p2 <- prod1Rang2 + 1
                p3 <- prod2Rang1 + 1
                p4 <- prod2Rang2 - 1
            }
            else {
                p1 <- prod1Rang1 + 1
                p2 <- prod1Rang2 - 1
                p3 <- prod2Rang1 - 1
                p4 <- prod2Rang2 + 1
            }
            ecap <- (p1 - unifPR)^2 + (p2 - unifPR)^2 + (p3 -
                unifPR)^2 + (p4 - unifPR)^2
            return(ecav - ecap)
        }
        nbEssPlan <- nbPanelist * nbRang
        weightSucc <- 1 - weightRang
        nbPertub <- 10
        matJROpt <- NA
        eamPondMin <- 1e+10
        matCond <- condenseRang(nbPanelist, nbProd, nbRang, nbPanelistImp,
            matJugeImp, plan)
        unifPP <- (nbPanelist * (nbRang - 1))/(nbProd * (nbProd -
            1))
        unifPR <- nbEssPlan/(nbProd * nbRang)
        for (iPlan in 1:nbDesign) {
            iPertu <- 1
            eamPondMinPertu <- 999999
            matJRPertu <- NA
            matPPPertu <- NA
            matPRPertu <- NA
            iStar <- -9
            prod1Max <- -1
            prod2Max <- -1
            rang1Max <- -1
            rang2Max <- -1
            prod1MaxAv <- -1
            prod2MaxAv <- -1
            rang1MaxAv <- -1
            rang2MaxAv <- -1
            iStarAv <- -1
            prod1MaxAnt <- -1
            prod2MaxAnt <- -1
            rang1MaxAnt <- -1
            rang2MaxAnt <- -1
            iStarAnt <- -1
            matJR <- genJRAlea(matCond, nbPanelist, nbRang, nbPanelistImp)
            matPP <- matrix(0, nbProd, nbProd)
            matPR <- matrix(0, nbProd, nbRang)
            for (i in 1:nbPanelist) {
                for (j in 2:nbRang) matPP[matJR[i, j - 1], matJR[i,
                  j]] <- matPP[matJR[i, j - 1], matJR[i, j]] +
                  1
                for (r in 1:nbRang) matPR[matJR[i, r], r] <- matPR[matJR[i,
                  r], r] + 1
            }
            it <- 0
            while (TRUE) {
                critMax <- -1000
                if (it >= 3) {
                  prod1MaxAnt <- prod1MaxAv
                  prod2MaxAnt <- prod2MaxAv
                  rang1MaxAnt <- rang1MaxAv
                  rang2MaxAnt <- rang2MaxAv
                  iStarAnt <- iStarAv
                }
                if (it >= 2) {
                  prod1MaxAv <- prod1Max
                  prod2MaxAv <- prod2Max
                  rang1MaxAv <- rang1Max
                  rang2MaxAv <- rang2Max
                  iStarAv <- iStar
                }
                valOpti <- valOptiPlan(nbPanelist, nbProd, nbRang,
                  nbEssPlan, unifPP, unifPR, weightRang, weightSucc,
                  matPR, matPP)
                if (valOpti < 1e-05)
                  break
                if (it >= 3 && iStar == iStarAnt && prod1Max ==
                  prod1MaxAnt && prod2Max == prod2MaxAnt && rang1Max ==
                  rang1MaxAnt && rang2Max == rang2MaxAnt) {
                  iStar <- iStar + 1
                  if (iPertu <= nbPertub) {
                    eamPR <- eam(matPR, unifPR)
                    eamPP <- eam(matPP, unifPP)
                    eamPondPertu <- weightRang * eamPR + weightSucc *
                      eamPP
                    if (eamPondPertu < eamPondMinPertu) {
                      matPRPertu <- matPR
                      matPPPertu <- matPP
                      matJRPertu <- matJR
                      eamPondMinPertu <- eamPondPertu
                    }
                    if (iStar <= nbPanelist) {
                      tmp <- matJR[iStar, 1]
                      for (j in 1:(nbRang - 1)) matJR[iStar,
                        j] <- matJR[iStar, j + 1]
                      matJR[iStar, nbRang] <- tmp
                    }
                    matPP[, ] <- 0
                    for (i in 1:nbPanelist) {
                      for (j in 2:nbRang) matPP[matJR[i, j -
                        1], matJR[i, j]] <- matPP[matJR[i, j -
                        1], matJR[i, j]] + 1
                    }
                    matPR[, ] <- 0
                    for (i in 1:nbPanelist) {
                      for (r in 1:nbRang) matPR[matJR[i, r],
                        r] <- matPR[matJR[i, r], r] + 1
                    }
                    iPertu <- iPertu + 1
                    next
                  }
                  else {
                    matPR <- matPRPertu
                    matJR <- matJRPertu
                    matPP <- matPPPertu
                    break
                  }
                }
                for (jugeActuel in (nbPanelistImp + 1):nbPanelist) {
                  for (rang1 in 1:(nbRang - 1)) {
                    for (rang2 in (rang1 + 1):nbRang) {
                      if (matJR[jugeActuel, rang1] < matJR[jugeActuel,
                        rang2]) {
                        prod1 <- matJR[jugeActuel, rang1]
                        prod2 <- matJR[jugeActuel, rang2]
                      }
                      else {
                        prod1 <- matJR[jugeActuel, rang2]
                        prod2 <- matJR[jugeActuel, rang1]
                      }
                      critPR <- critQualPR(prod1, rang1, prod2,
                        rang2, matPR, unifPR, matJR, jugeActuel)
                      critPP <- critQualPP(rang1, rang2, nbRang,
                        matPP, unifPP, matJR, jugeActuel)
                      critGlob <- weightRang * critPR + weightSucc *
                        critPP
                      if (critGlob > critMax) {
                        critMax <- critGlob
                        iStar <- jugeActuel
                        prod1Max <- prod1
                        prod2Max <- prod2
                        rang1Max <- rang1
                        rang2Max <- rang2
                      }
                    }
                  }
                }
                matPP <- majPP(matPP, iStar, matJR, nbRang, rang1Max,
                  rang2Max)
                matPR <- majPR(matPR, rang1Max, rang2Max, iStar,
                  matJR, prod1Max, prod2Max)
                matJR <- majJR(matJR, rang1Max, rang2Max, iStar)
                it <- it + 1
            }
            eamPR <- eam(matPR, unifPR)
            eamPP <- eam(matPP, unifPP)
            eamPond <- weightRang * eamPR + weightSucc * eamPP
            if (eamPond < eamPondMin) {
                matJROpt <- matJR
                eamPondMin <- eamPond
            }
        }
        return(matJROpt)
    }
    ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
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
    fedorov <- function(nbMod1, nbMod2, nbEssPlan, matCond, matEffet,
        tabEssais, nbEssImp = 0, nbDesign = 10) {
        planValide <- function(plan, npj) {
            test <- unique(rowSums(plan))
            if (length(test) != 1 | test[1] != npj)
                return(FALSE)
            else return(TRUE)
        }
        genJPAlea <- function(nbEssaiImp, indEssImp, nbEssaiPlan,
            nbEssaiPlanComp) {
            res <- vector(mode = "integer", length = nbEssaiPlanComp)
            dejaTire <- vector(length = nbEssaiPlanComp)
            if (nbEssaiImp > 0) {
                res[1:nbEssaiImp] <- indEssImp[1:nbEssaiImp]
                dejaTire[res[1:nbEssaiImp]] <- TRUE
            }
            nbTirage <- nbEssaiPlan - nbEssaiImp
            i <- 0
            while (i < nbTirage) {
                tirage <- as.integer(1 + runif(1) * nbEssaiPlanComp)
                if (dejaTire[tirage])
                  next
                res[nbEssaiImp + 1 + i] <- tirage
                dejaTire[tirage] <- TRUE
                i <- i + 1
            }
            return(res)
        }
        tabEssaisCpy <- tabEssais
        nbEssTotal <- nbMod1 * nbMod2
        modTotalEff <- (nbMod1 + nbMod2 - 1)
        planOpt <- matrix()
        detOpt <- -9e+10
        iPlan <- 0
        nbMod1ParMod2 <- nbEssPlan/nbMod1
        if (nbEssImp + 1 >= nbEssPlan) {
            nbPres <- vector(mode = "integer", length = nbEssTotal)
            nbPres[tabEssaisCpy[1:nbEssPlan]] <- nbPres[tabEssaisCpy[1:nbEssPlan]] +
                1
            res <- matrix(0, nrow = nbMod1, ncol = nbMod2)
            for (p in 1:nbEssTotal) res[matCond[p, 1], matCond[p,
                2]] <- nbPres[p]
            return(res)
        }
        while (iPlan < nbDesign) {
            matEffetSub <- matEffet[tabEssaisCpy[1:nbEssPlan],
                ]
            PX <- t(matEffetSub) %*% matEffetSub
            detPlan <- det(PX, logarithme = FALSE)
            if (detPlan < 1e-05) {
                tabEssaisCpy <- genJPAlea(nbEssaiImp = nbEssImp,
                  indEssImp = tabEssaisCpy, nbEssaiPlan = nbEssPlan,
                  nbEssaiPlanComp = nbEssTotal)
                next
            }
            DISP <- ginv(PX)
            nbPres <- vector(mode = "integer", length = nbEssTotal)
            nbPres[tabEssaisCpy[1:nbEssPlan]] <- nbPres[tabEssaisCpy[1:nbEssPlan]] +
                1
            while (TRUE) {
                jmax <- 0
                kmax <- 0
                deltaMax <- 0
                tabVar <- diag((matEffet %*% DISP) %*% t(matEffet))
                if (nbEssImp + 1 >= nbEssPlan)
                  break
                for (k in c((nbEssImp + 1):nbEssPlan)) {
                  essaiK <- tabEssaisCpy[k]
                  covar1 <- (matEffet %*% t(DISP))[essaiK, ]
                  for (j in 1:nbEssTotal) {
                    if (nbPres[j] < 0.5 && tabVar[j] > tabVar[essaiK]) {
                      covkj <- sum(covar1 %*% (t(matEffet))[,
                        j])
                      delta <- tabVar[j] * (1 - tabVar[essaiK]) -
                        tabVar[essaiK] + covkj^2
                      if (delta > deltaMax) {
                        deltaMax <- delta
                        jmax <- j
                        kmax <- k
                      }
                    }
                  }
                }
                if (deltaMax < 1e-15)
                  break
                Z <- (DISP %*% t(matEffet))[, jmax]
                ZZ <- matrix(0, modTotalEff, modTotalEff)
                for (i2 in 1:modTotalEff) {
                  for (i1 in 1:modTotalEff) ZZ[i1, i2] <- DISP[i1,
                    i2] - ((Z[i1] * Z[i2])/(1 + tabVar[jmax]))
                }
                matEffetSub <- matEffet[tabEssaisCpy[kmax], ]
                vart <- sum((matEffetSub %*% t(ZZ)) %*% matEffetSub)
                Z <- t(ZZ %*% t(matEffet)[, tabEssaisCpy[kmax]])
                for (i2 in 1:modTotalEff) {
                  for (i1 in 1:modTotalEff) DISP[i1, i2] <- ZZ[i1,
                    i2] + (Z[i1] * Z[i2])/(1 - vart)
                }
                nbPres[jmax] <- nbPres[jmax] + 1
                nbPres[tabEssaisCpy[kmax]] <- nbPres[tabEssaisCpy[kmax]] -
                  1
                tabEssaisCpy[kmax] <- jmax
                detPlan <- detPlan * (1 + deltaMax)
            }
            res <- matrix(data = 0, nrow = nbMod1, ncol = nbMod2)
            for (p in 1:nbEssTotal) res[matCond[p, 1], matCond[p,
                2]] <- nbPres[p]
            if (!planValide(res, nbMod1ParMod2)) {
                iPlan <- iPlan + 1
                if (iPlan >= nbDesign) {
                  print(res)
                  stop("\nNo solution: \nTry again or increase the number of judges or the number of products per judge.\n")
                }
                else next
            }
            if (detPlan > detOpt) {
                detOpt <- detPlan
                planOpt <- res
            }
            tabEssaisCpy <- genJPAlea(nbEssaiImp = nbEssImp,
                indEssImp = tabEssais, nbEssaiPlan = nbEssPlan,
                nbEssaiPlanComp = nbEssTotal)
            iPlan <- iPlan + 1
        }
        return(planOpt)
    }
    if (nbPanelist < nbPanelistMin || nbProdByPanelist > nbProd)
        stop("CPDO : problem in the parameters!")
    set.seed(graine)
    res <- list(design = list(), rank = list())
    nbEssJuge <- 0
    matEssJuge <- NA
    nbEssImp <- 0
    if (!is.na(matEssImp))
        nbEssImp <- nrow(matEssImp)
    for (nbPanelistActuel in (nbPanelistMin:nbPanelist)) {
        nep <- nbPanelistActuel * nbProdByPanelist
        nep.complet <- nbPanelistActuel * nbProd
        if (nbPanelistActuel + nbProd - 1 > nep)
            stop("Matrix not inversible: not sufficiently experiments")
        matCond <- condenseProd(nbPanelist = nbPanelistActuel,
            nbProd = nbProd, nbEssai = nep.complet)
        matCandidate <- matEffet(nbPanelist = nbPanelistActuel,
            nbProd = nbProd, matCond = matCond)
        if (nbEssImp > 0)
            tabEssais <- findEssImp(matEssImp, matCond)
        else tabEssais <- NULL
        if (nbEssImp < nep)
            tabEssais <- genJPAlea(nbEssaiImp = nbEssImp, indEssImp = tabEssais,
                nbEssaiPlan = nep, nbEssaiPlanComp = nep.complet)
        planOptProd <- fedorov(nbPanelistActuel, nbProd, nep,
            matCond, matCandidate, tabEssais, nbEssImp, nbDesign = nbDesignProd)
        rownames(planOptProd) = paste("Panelist", 1:nrow(planOptProd))
        colnames(planOptProd) = paste("Product", 1:ncol(planOptProd))
        lstEssais <- calcEssaisImp(planOptProd)
        nbEssImp <- lstEssais$nb
        matEssImp <- lstEssais$valeurs
        planOptOrdre <- NA
        if (ordre) {
            planOptOrdre <- rang(nbPanelistActuel, nbProd, nbProdByPanelist,
                planOptProd, weight, nbPanelistImp = nbEssJuge,
                matJugeImp = matEssJuge, nbDesignOrdre)
            nbEssJuge <- nrow(planOptOrdre)
            matEssJuge <- planOptOrdre
            rownames(planOptOrdre) = paste("Panelist", 1:nrow(planOptOrdre))
            colnames(planOptOrdre) = paste("Rank", 1:ncol(planOptOrdre))
        }
    }
    res$design <- planOptProd
    res$rank <- planOptOrdre
    class(res) <- list("planDO", "list")
    return(res)
}
