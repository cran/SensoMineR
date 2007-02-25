# Construit des plans de dégustation optimaux à l'aide de l'algorithme
# d'échange de Fedorov, en tenant compte de quatre effets : juge, produit,
# rang et succession. Cette fonction permet également de construire un plan
# à partir d'un nombre "incertain" de juges.
# @nbPanelist        : Nombre total de juges.
# @nbPanelistMin     : Nombre de juges minimum.
# @nbProd        : Nombre total de produits.
# @nbProdByPanelist : Nombre de produits à attribuer à chaque juge.
# @ordre         : Prise en compte de l'ordre (TRUE ou FALSE).
# @weight         : Importance relative rang/succesion (0.0 à 1.0).
# @gaine         : Graine du générateur aléatoire.
# @nbDesignProd    : Nombre de plans à construire pour l'affectation des produits.
# @nbDesignOrdre   : Nombre de plans à construire pour l'affectation de l'ordre
#                  de présentation des produits.
# @nbEssaiImp    : Nombre d'essais imposés.
# @matEssaiImp   : Tableau contenant les essais imposés.
# return         : Une liste où l'élément $plan contient les plans optimaux
#                  et l'élément $ordre contient les ordres de présentation
#                  des produits.
# ========================================================================== #

optimaldesign<-function(nbPanelist, nbProd, nbProdByPanelist, nbPanelistMin = nbPanelist,  ordre = TRUE, weight = 0.5, graine = Sys.time(), nbDesignProd = 10, nbDesignOrdre = 50, matEssImp = NA )
{

  # Fonctions internes
  # La fonction condenseProd crée une matrice condensée Juge*produit, Toutes les
  # combinaison possible entre les facteurs juge et produit sont créées et ordonnées
  # par juge puis par produit.
  # L'argument nbEssai permet de ne retourner que les nbEssai premières lignes
  condenseProd<-function (nbPanelist, nbProd, nbEssai=NULL)  {
    colJuge<-gl(nbPanelist, nbProd)
    colProd<-as.factor(seq(nbProd))
    if(is.null(nbEssai) | nbEssai>(nbPanelist*nbProd)) nbEssai<-nbPanelist*nbProd
    res<-cbind.data.frame(colJuge,colProd)
    # On retourne les nbEssai premiers essais demandés
    return(res[1:nbEssai, ])
  }

  # Création à partir d'un plan d'expérience présenté sous forme de sa matrice
  # condensée, de la matrice des effets associé ( = matrice binaire à nbre d'essais
  # lignes et à nbre total de modalité +1 colonnes
  matEffet<- function(nbPanelist, nbProd, matCond) {
    res<-model.matrix(~colProd+colJuge, data=matCond, contrasts=list(colProd="contr.sum", colJuge="contr.sum"))
    # Réorganisation de la matrice : juge1...jugeN ; prod1...prodP ; CSTE
    partieProd<-2:nbProd
    partieJuge<-c((nbProd+1) : (nbProd+nbPanelist-1))
    res<-res[ , c(partieJuge, partieProd, 1)]
    return (res)
  }

  # nbEssaiImp : Nombre d'essais imposés
  # nbEssaiPlan : Nombre d'essais que l'on souhaite conserver
  # nbEssai : Nombre d'essais initiale
  genJPAlea<-function(nbEssaiImp, indEssImp, nbEssaiPlan, nbEssaiPlanComp) {
    res<-vector(mode="integer", length=nbEssaiPlanComp)
    dejaTire<-vector(length=nbEssaiPlanComp)

    if (nbEssaiImp>0) {
      # Recopie des anciennes données
      res[1:nbEssaiImp]<-indEssImp[1:nbEssaiImp]
      # On sauvegarde les numéros d'essais imposés d'indice 1 à nbEssaiImp
      dejaTire[res[1:nbEssaiImp]]<-TRUE
    }

    # Tirage aléatoire ...
    nbTirage<-nbEssaiPlan-nbEssaiImp
    i<-0
    while (i < nbTirage) {
      tirage<-as.integer(1+runif(1)*nbEssaiPlanComp)
      if (dejaTire[tirage]) next
      res[nbEssaiImp+1+i]<-tirage
      dejaTire[tirage]<-TRUE
      i<-i+1
    }
    return (res)
  }

  # Crée une matrice contenant les essais imposés (juge i * prod j) à partir
  # d'un plan optimal.
  # @plan  : Plan d'expérience.
  # return : Liste contenant le nombre d'essais imposés et la matrice.
  calcEssaisImp<-function(plan) {
    nbLig<-nrow(plan)
    nbCol<-ncol(plan)
    matEssImp<-matrix(nrow=nbLig*nbCol, ncol=2)
    k <- 0
    for (i in 1:nbLig) {
      for (j in 1:nbCol) {
        if (plan[i,j] == 1) {
          k<-k+1
          matEssImp[k,1]<-i
          matEssImp[k,2]<-j
        }
      }
    }

    return(list(nb=k, valeurs=matEssImp[1:k, ]))
  }

  # Retourne un vecteur correspondant aux N° de lignes des essais imposés dans
  # la matrice condensée du plan complet
  findEssImp<-function(matEssImp, matCond) {
    nbEssImp = nrow(matEssImp)
    res<-vector(mode="integer", length=nbEssImp)
    for (ei in 1:nbEssImp) res[ei]<-which(matCond[,1]==matEssImp[ei,1] & matCond[,2]==matEssImp[ei,2])
    return (res)
  }

######### Fonctions pour calculer les rangs
# optimisation de l'ordre de presentation d'un plan (Juge x Produit) 
# sur les critère ordre et succession. On peut également spécifier le poids
# des effets produit-produit et produit-rang que l'on souhaite imposer.
# @nbPanelist     : Nombre total de juges.
# @nbProd     : Nombre total de produits.
# @nbRang     : Nombre de produis à attribuer à chaque juge.
# @plan       : Plan JugeXProduit optimal.
# @weightRang  : Importance relative rang/succesion (0.0 à 1.0).
# @nbPanelistImp  : Nombre de juges imposé.
# @matJugeImp : Matrice contenant les juges imposés.
# @nbDesign     : Nombre de plans à construire.
# @nbIter     : Nombre maximum d'itération de la procédure d'échange.
# return      : Un ordre de présentation des produits.
# ========================================================================== #

rang<-function(nbPanelist, nbProd, nbRang, plan, weightRang, nbPanelistImp, matJugeImp, nbDesign)
{
  # fonctions internes
  # Crée une matrice Juge x Rang aléatoirement. Il est possible d'imposer des
  # juges et des rangs.
  # @matCond   : Matrice condensée.
  # @nbPanelist    : Nombre de juges.
  # @nbRang    : Nombre de rang.
  # @nbPanelistImp : Nombre de juges imposé.
  # return     : Une matrice Juge x Rang aléatoirement générer.
  genJRAlea<-function(matCond, nbPanelist, nbRang, nbPanelistImp)
  {
    # Matrice de travail
    matAlea<-matCond

    # Génération aléatoire
    for (i in (nbPanelistImp+1):nbPanelist) {
      tot<-nbRang
      for (r in 1:(nbRang-1)) {
        jal<-as.integer(1+runif(1)*tot)
        matAlea[i,r]<-matCond[i,jal]
        if (jal<=tot-1) {
          for (k in jal:(tot-1)) matCond[i,k]<-matCond[i,k+1]
        }
        tot<-tot-1
      }
      matAlea[i,nbRang]<-matCond[i,1]
    }
    return (matAlea)
  }


  # Crée une matrice condensée contenant pour chaque juge, l'indice du produit
  # dégusté dans le plan PxJ issu de Fedorov.
  # @nbPanelist     : Nombre de juges.
  # @nbProd     : Nombre de produits.
  # @nbRang     : Nombre de rang.
  # @nbPanelistImp  : Nombre de juges imposés.
  # @matJugeImp : Matrice contenant les juges imposés.
  # @plan       : Plan Juge x Produit issus de Fedorov.
  # @return     : Une matrice condensée.
  condenseRang <- function (nbPanelist, nbProd, nbRang, nbPanelistImp, matJugeImp, plan)
  {
    #vérification des paramètres
    if(nbPanelistImp!=0) {
      if (nbPanelistImp!=nrow(matJugeImp)) stop("renseignement des juges imposés incohérent")
    }
    jugeProCond<-matrix(0, nbPanelist, nbRang)

    if(nbPanelistImp==0) {
      for (i in 1:nbPanelist) {
        k<-1
        for (j in 1:nbProd) {
          if (plan[i,j]==1)  {
            jugeProCond[i,k]<-j
            k<-k+1
          }
        }
      }
    }
    else {
      jugeProCond[1:nbPanelistImp,]<-matJugeImp
      for(i in ((nbPanelistImp+1):nbPanelist)) {
        k<-1
        for (j in 1:nbProd) {
          if (plan[i,j]==1) {
            jugeProCond[i,k]<-j
            k<-k+1
          }
        }
      }
    }
    return (jugeProCond)
  }


  # Calcule l'Ecart Absolu Moyen (EAM) d'une matrice par rapport à son effectif
  # théorique.
  # @mat   : Matrice.
  # @effTh : Effectif théorique.
  # return : L'écart absolu moyen de 'mat'.
  eam<-function(mat, effTh)
  {
    ecartUnif<-sum(abs(effTh-mat))
    return(ecartUnif/(nrow(mat)*ncol(mat)))
  }

  # Echange dans la matrice Juge x Rang le produit du rang1 avec celui du rang2
  # pour le juge i.
  # @matJR    : Matrice JxR sur laquelle on applique l'échange.
  # @rang1Max : Rang du produit1 à échanger.
  # @rang2Max : Rang du produit2 à échanger.
  # @i        : Juge sur lequel on travaille.
  # return    : La matrice JxR mise à jour avec les échanges.
  majJR<-function(matJR, rang1Max, rang2Max, i)
  {
    tmp<-matJR[i,rang1Max]
    matJR[i,rang1Max]<-matJR[i,rang2Max]
    matJR[i,rang2Max]<-tmp

    return (matJR)
  }

  # Met à jour les successions dans la matrice PxP après une échange du
  # produit de rang1 et du produit de rang2 pour un juge donné.
  # @matPP    : Matrice PxP sur laquelle on applique l'échange.
  # @i        : Juge sur lequel on travaille.
  # @matJR    : Matrice Juge x Rang.
  # @nbRang   : Nombre de rang.
  # @rang1Max : Rang du produit1 à échanger.
  # @rang2Max : Rang du produit2 à échanger.
  # return    : La matrice PxP mise à jour avec les échanges.
  majPP<-function(matPP, i, matJR, nbRang, rang1Max, rang2Max)
  {

    # MAJ des (-1)
    if (rang1Max != 1) matPP[matJR[i,rang1Max-1],matJR[i,rang1Max]]<-matPP[matJR[i,rang1Max-1],matJR[i,rang1Max]] -1
    if (rang1Max+1 == rang2Max) matPP[matJR[i,rang1Max],matJR[i,rang2Max]]<-matPP[matJR[i,rang1Max],matJR[i,rang2Max]] -1
    else {
      matPP[matJR[i,rang1Max],matJR[i,rang1Max+1]]<-matPP[matJR[i,rang1Max],matJR[i,rang1Max+1]] -1
      matPP[matJR[i,rang2Max-1],matJR[i,rang2Max]]<-matPP[matJR[i,rang2Max-1],matJR[i,rang2Max]] -1
    }

    if (rang2Max != nbRang) matPP[matJR[i,rang2Max],matJR[i,rang2Max+1]]<-matPP[matJR[i,rang2Max],matJR[i,rang2Max+1]] -1

    # MAJ des (+1)
    if (rang2Max != nbRang) matPP[matJR[i,rang1Max],matJR[i,rang2Max+1]]<-matPP[matJR[i,rang1Max],matJR[i,rang2Max+1]] +1
    if (rang1Max != 1) matPP[matJR[i,rang1Max-1],matJR[i,rang2Max]]<-matPP[matJR[i,rang1Max-1],matJR[i,rang2Max]] +1
    if (rang2Max != rang1Max+1) {
      matPP[matJR[i,rang2Max],matJR[i,rang1Max+1]]<-matPP[matJR[i,rang2Max],matJR[i,rang1Max+1]] +1
      matPP[matJR[i,rang2Max-1],matJR[i,rang1Max]]<-matPP[matJR[i,rang2Max-1],matJR[i,rang1Max]] +1
    }
    else matPP[matJR[i,rang2Max],matJR[i,rang1Max]]<-matPP[matJR[i,rang2Max],matJR[i,rang1Max]] +1
    return (matPP)
  }

  # @matPR    : Matrice PxP sur laquelle on applique l'échange.
  # @rang1Max : Rang du produit1 à échanger.
  # @rang2Max : Rang du produit2 à échanger.
  # @i        : Juge sur lequel on travaille.
  # @matJR    : Matrice Juge x Rang.
  # @prod1Max : Produit1 sur lequel, on a modifié son rang.
  # @prod2Max : Produit2 sur lequel, on a modifié son rang.
  # return    : La matrice PxR mise à jour avec les échanges.
  majPR<-function( matPR, rang1Max, rang2Max, i, matJR, prod1Max, prod2Max )
  {
    if (matJR[i,rang1Max]<matJR[i,rang2Max]) {
      matPR[prod1Max,rang1Max]<-matPR[prod1Max,rang1Max] -1
      matPR[prod1Max,rang2Max]<-matPR[prod1Max,rang2Max] +1
      matPR[prod2Max,rang1Max]<-matPR[prod2Max,rang1Max] +1
      matPR[prod2Max,rang2Max]<-matPR[prod2Max,rang2Max] -1
    }
    else {
      matPR[prod1Max,rang1Max]<-matPR[prod1Max,rang1Max] +1
      matPR[prod1Max,rang2Max]<-matPR[prod1Max,rang2Max] -1
      matPR[prod2Max,rang1Max]<-matPR[prod2Max,rang1Max] -1
      matPR[prod2Max,rang2Max]<-matPR[prod2Max,rang2Max] +1
    }

    return (matPR)
  }

  # Calcule le nombre théorique moyen d'apparition d'un produit pour chacun des
  # rangs et du nombre de couples de produits ordonnés (produit précédent -
  # produit suivant).
  # @nbPanelist    : Nombre total de juges.
  # @nbProd    : Nombre total de produits.
  # @nbRang    : Nombre de produis à attribuer à chaque juge.
  # @nbEssPlan : Nombre d'essais du plan.
  # @unifPP    : Effectif théorique pour PxP.
  # @unifPR    : Effectif théorique pour PxR.
  # @weightRang : Poids de la maitrise de l'effet rang.
  # @weightSucc : Poids de la maitrise de l'effet succession.
  # @matPR     : Matrice Produit x Rang.
  # @matPP     : Matrice Produit x Produit.
  # return     : Le delta entre la théorie et la réalité.
  valOptiPlan<-function(nbPanelist, nbProd, nbRang, nbEssPlan, unifPP, unifPR, weightRang, weightSucc, matPR, matPP )
  {
    # Optimalité Rang
    optiRang1<-sum(abs((nbEssPlan/nbRang)-colSums(matPR)))
    optiRang2<-sum(as.integer(abs(unifPR-matPR)))
    optiRang<-optiRang1+optiRang2


    # Optimalité Succ
    totslig<-vector(mode="integer", length=nbProd)
    for (i in 1:nbProd) {
      for (j in 1:nbProd) {
        if(i!=j) totslig[i]<-totslig[i]+matPP[i,j]
      }
    }

    optiSucc1<-sum(as.integer (abs(((nbEssPlan-nbPanelist)/nbProd) - totslig)))
    optiSucc2<-(abs (unifPP-matPP))
    diag(optiSucc2)<-0
    optiSucc2<-sum(as.integer(optiSucc2))
    optiSucc<-optiSucc1+optiSucc2

    return(weightRang*optiRang+weightSucc*optiSucc)
  }

  # Calcule la qualité du plan PxP après une modification de la succession de
  # deux produits donnés.
  # @rang1  : Rang du produit1.
  # @rang2  : Rang du produit2.
  # @nbRang : Nombre de rang (= nombre de produit/juge).
  # @matPP  : Matrice produit/produit (PP).
  # @unifPP : Uniformité de la matrice PxP.
  # @matJR  : Matrice Juge-Rang
  # @juge   : Juge sur lequel on travaille dans la matrice JxR.
  # return  : La qualité de l'échange.
  critQualPP<-function(rang1, rang2, nbRang, matPP, unifPP, matJR, juge)
  {
    # Distance avant modification
    distuv<-0
    # Distance après modification
    distup<-0

    # ----------------------------------------------------
    # Calcul de la distance à l'uniformité avant l'échange
    # ----------------------------------------------------
    # Calcul des termes en (-)
    if (rang1 != 1) distuv <- distuv + (matPP[matJR[juge,rang1-1],matJR[juge,rang1]] - unifPP) ^2
    if (rang1+1 == rang2) {
      distuv <- distuv + (matPP[matJR[juge,rang1],matJR[juge,rang2]] - unifPP) ^2
      distuv <- distuv + (unifPP - matPP[matJR[juge,rang2],matJR[juge,rang1]]) ^2
    }
    else {
      distuv <- distuv + (matPP[matJR[juge,rang1],matJR[juge,rang1+1]] - unifPP) ^2
      distuv <- distuv + (matPP[matJR[juge,rang2-1],matJR[juge,rang2]] - unifPP) ^2
    }
    if (rang2 != nbRang) distuv <- distuv + (matPP[matJR[juge,rang2],matJR[juge,rang2+1]] - unifPP) ^2

    # Calcul des termes en (+)
    if (rang1+1 != rang2) distuv <- distuv + (unifPP - matPP[matJR[juge,rang2-1],matJR[juge,rang1]]) ^2
    if (rang2 != nbRang) distuv <- distuv + (unifPP - matPP[matJR[juge,rang1],matJR[juge,rang2+1]]) ^2
    if (rang1 != 1) distuv <- distuv + (unifPP - matPP[matJR[juge,rang1-1],matJR[juge,rang2]]) ^2
    if (rang2 != rang1+1) distuv <- distuv + (unifPP - matPP[matJR[juge,rang2],matJR[juge,rang1+1]]) ^2


    # ----------------------------------------------------
    # Calcul de la distance à l'uniformité après l'échange
    # ----------------------------------------------------
    # Calcul des termes en (-)
    if (rang1 != 1) distup <- distup + (unifPP - matPP[matJR[juge,rang1-1],matJR[juge,rang1]] +1) ^2
    if (rang1+1 == rang2) {
      distup <- distup + (unifPP - matPP[matJR[juge,rang1],matJR[juge,rang2]]+1) ^2
      distup <- distup + (unifPP - matPP[matJR[juge,rang2],matJR[juge,rang1]]-1) ^2
    }
    else {
      distup <- distup + (unifPP - matPP[matJR[juge,rang1],matJR[juge,rang1+1]] +1) ^2
      distup <- distup + (unifPP - matPP[matJR[juge,rang2-1],matJR[juge,rang2]] +1) ^2
    }
    if (rang2 != nbRang) distup <- distup + (unifPP - matPP[matJR[juge,rang2],matJR[juge,rang2+1]]+1) ^2

    # Calcul des termes en (+)
    if (rang1+1 != rang2) distup <- distup + (unifPP - matPP[matJR[juge,rang2-1],matJR[juge,rang1]]-1) ^2
    if (rang2 != nbRang) distup <- distup + (unifPP - matPP[matJR[juge,rang1],matJR[juge,rang2+1]]-1) ^2
    if (rang1 != 1) distup <- distup + (unifPP - matPP[matJR[juge,rang1-1],matJR[juge,rang2]]-1) ^2
    if (rang2 != rang1+1) distup <- distup + (unifPP - matPP[matJR[juge,rang2],matJR[juge,rang1+1]]-1) ^2
    return (distuv - distup)
  }

  # Calcule la qualité du plan PxR après une modification du rang de deux
  # produits donnés.
  # @prod1  : Produit 1.
  # @prod2  : Produit 2.
  # @rang1  : Rang du produit1.
  # @rang2  : Rang du produit2.
  # @nbRang : Nombre de rang (= nombre de produit/juge).
  # @matPR  : Matrice produit/rang (PxR).
  # @unifPR : Uniformité de la matrice PxR.
  # @matJR  : Matrice Juge-Rang
  # @juge   : Juge sur lequel on travaille dans la matrice JxR.
  # return  : La qualité de l'échange.
  critQualPR<-function(prod1, rang1, prod2, rang2, matPR, unifPR, matJR, juge)
  {
    # Nombre de fois qu'un produit est dégusté à un certain rang
    prod1Rang1<-matPR[prod1,rang1]
    prod1Rang2<-matPR[prod1,rang2]
    prod2Rang1<-matPR[prod2,rang1]
    prod2Rang2<-matPR[prod2,rang2]

    # Ecart à son effectif uniforme avant l'échange de test
    ecav<-(prod1Rang1 - unifPR)^2 + (prod1Rang2 - unifPR)^2 + (prod2Rang1 - unifPR)^2 + (prod2Rang2 - unifPR)^2

    # Echange de test
    if (matJR[juge,rang1] < matJR[juge,rang2]) {
      p1<-prod1Rang1 - 1
      p2<-prod1Rang2 + 1
      p3<-prod2Rang1 + 1
      p4<-prod2Rang2 - 1
    }
    else {
      p1<-prod1Rang1 + 1
      p2<-prod1Rang2 - 1
      p3<-prod2Rang1 - 1
      p4<-prod2Rang2 + 1
    }

    # Ecart à son effectif uniforme après l'échange de test
    ecap<-(p1 - unifPR)^2 + (p2 - unifPR)^2 + (p3 - unifPR)^2 + (p4 - unifPR)^2

    return (ecav - ecap)
  }

# ========================================================================== #
# fin des fonctions internes
# ========================================================================== #



  nbEssPlan<-nbPanelist*nbRang  # Nombre d'essais total du plan
  weightSucc<-1-weightRang  # Poids succession des produits
  nbPertub<-10  # Nombre de perturbation du plan
  matJROpt<-NA  # definition de l'objet matrice (Juge x rang) à optimiser
  eamPondMin<-1E10  # EAM pondéré du meilleur plan mémorisé


  # Matrice condensée (JxR) du plan optimal issus de Fedorov
  matCond<-condenseRang(nbPanelist, nbProd, nbRang, nbPanelistImp, matJugeImp, plan) 
  
  unifPP<-(nbPanelist*(nbRang-1))/(nbProd*(nbProd-1))  # Nombre de fois qu'un couple de produits peut apparaitre
  unifPR<-nbEssPlan/(nbProd*nbRang)  # Nombre de fois qu'un couple produit-rang peut apparaitre


  # Construction de 'nbPlan' plans de rang optimaux
  for (iPlan in 1:nbDesign) {
    iPertu<-1      # Compteur du nombre de perturbation effectuée
    eamPondMinPertu<-999999      # Meilleur EAM pondéré après perturbation
    matJRPertu<-NA    # Meilleure matrice JxR perturbée
    matPPPertu<-NA    # Meilleure matrice PxP perturbée
    matPRPertu<-NA    # Meilleure matrice PxR perturbée
    iStar<- -9    # Indice du juge qui a subit un échange
    

    # Contient les meilleurs paramètres de l'échange actuel
    prod1Max<- -1
    prod2Max<- -1
    rang1Max<- -1
    rang2Max<- -1

    # Contient les meilleurs paramètres du dernier échange
    prod1MaxAv<- -1
    prod2MaxAv<- -1
    rang1MaxAv<- -1
    rang2MaxAv<- -1
    iStarAv<- -1

    # Contient les meilleurs paramètres de l'avant dernier échange
    prod1MaxAnt<- -1
    prod2MaxAnt<- -1
    rang1MaxAnt<- -1
    rang2MaxAnt<- -1
    iStarAnt<- -1

    # Construction aléatoire de la matrice JxR
    matJR<-genJRAlea(matCond, nbPanelist, nbRang, nbPanelistImp)
    # Construction de la matrice PxP qui en découle
    # Construction de la matrice PxR qui en découle
    matPP<-matrix(0, nbProd, nbProd)
    matPR<-matrix(0, nbProd, nbRang)
    for (i in 1:nbPanelist)  {
      for (j in 2:nbRang)  matPP[matJR[i,j-1],matJR[i,j]] <- matPP[matJR[i,j-1],matJR[i,j]] + 1
      for (r in 1:nbRang) matPR[matJR[i,r],r]<-matPR[matJR[i,r],r]+1
    }

    # Itération de la procédure d'échange
    it<-0
    while(TRUE) {
      critMax<- -1000

      # Sauvegarde des produits et rangs optimaux de l'avant dernier échange
      if (it>=3) {
        prod1MaxAnt<-prod1MaxAv
        prod2MaxAnt<-prod2MaxAv
        rang1MaxAnt<-rang1MaxAv
        rang2MaxAnt<-rang2MaxAv
        iStarAnt<-iStarAv
      }
      
      # Sauvegarde des produits et rangs optimaux du dernier échange
      if (it>=2) {
        prod1MaxAv<-prod1Max
        prod2MaxAv<-prod2Max
        rang1MaxAv<-rang1Max
        rang2MaxAv<-rang2Max
        iStarAv<-iStar
      }

      # Arret si la valeur optimale du plan est égale à 0
      valOpti<-valOptiPlan(nbPanelist, nbProd, nbRang, nbEssPlan, unifPP, unifPR, weightRang, weightSucc, matPR, matPP)
      if (valOpti<1E-5) break

      # Arret si les couples de produits sont identiques depuis les 2
      # derniers échanges
      if (it>=3 && iStar==iStarAnt && prod1Max==prod1MaxAnt && prod2Max==prod2MaxAnt && rang1Max==rang1MaxAnt && rang2Max==rang2MaxAnt) {
       # Permet d'éviter de satisfaire le test ci-dessus
       iStar<-iStar+1
        
       if (iPertu<=nbPertub) {
         # Calcul de l'EAM (Ecart Absolu Moyen) associé à PxR
         eamPR<-eam(matPR, unifPR)
         # Calcul de l'EAM associé à PxP
         eamPP<-eam(matPP, unifPP)
         
         # Calcul de l'EAM pondéré
         eamPondPertu<-weightRang*eamPR+weightSucc*eamPP
         if (eamPondPertu<eamPondMinPertu) {
           matPRPertu<-matPR
           matPPPertu<-matPP
           matJRPertu<-matJR
           eamPondMinPertu<-eamPondPertu
         }

         # Perturbation des rangs dans JxR
## modif Husson         
         if (iStar<=nbPanelist) {
           tmp<-matJR[iStar,1]
           for (j in 1:(nbRang-1)) matJR[iStar,j]<-matJR[iStar,j+1]
           matJR[iStar,nbRang]<-tmp
         }
#         if (iStar<=nbPanelist) matJR[iStar,]=matJR[iStar,c(nbRang,1:(nbRang-1))]
          
         # Mise à jour sur PxP
         matPP[,]<-0
         for (i in 1:nbPanelist) {
           for (j in 2:nbRang) matPP[matJR[i,j-1],matJR[i,j]]<-matPP[matJR[i,j-1],matJR[i,j]]+1
         }
          
         # Mise à jour de la matrice PxR
         matPR[,]<-0
         for (i in 1:nbPanelist) {
           for (r in 1:nbRang)  matPR[matJR[i,r],r]<-matPR[matJR[i,r],r] + 1
         }

         iPertu<-iPertu+1
         next    
       }
       
       # On recupère les meilleurs plans perturbés
       else {
         matPR<-matPRPertu
         matJR<-matJRPertu
         matPP<-matPPPertu
         break
       }
      }

      # On recherche le meilleur couple de produits à échanger
      for (jugeActuel in (nbPanelistImp + 1):nbPanelist) {
        for (rang1 in 1:(nbRang-1)) {
          for (rang2 in (rang1+1):nbRang) {
            # Sélection de 2 produits à échanger
            if (matJR[jugeActuel,rang1]<matJR[jugeActuel,rang2])  {
              prod1<-matJR[jugeActuel,rang1]
              prod2<-matJR[jugeActuel,rang2]
            }
            else {
              prod1<-matJR[jugeActuel,rang2]
              prod2<-matJR[jugeActuel,rang1]
            }
            
            # Calcul de la qualité du quadruplet associé dans PxR
            critPR<-critQualPR(prod1, rang1, prod2, rang2, matPR, unifPR, matJR, jugeActuel)

            # Calcul de la qualité du quadruplet associé dans PxP
            critPP<-critQualPP(rang1, rang2, nbRang, matPP, unifPP, matJR, jugeActuel)

            # Sauvegarde des meilleurs couples...
            critGlob<-weightRang*critPR+weightSucc*critPP
            if (critGlob>critMax) {
              critMax<-critGlob
              iStar<-jugeActuel
              prod1Max<-prod1
              prod2Max<-prod2
              rang1Max<-rang1
              rang2Max<-rang2
            }
          }
        }  
      }


#### Début  modif Husson         
##jugeActuel <- (nbPanelistImp + 1):nbPanelist
##rang1 <- rep(1:(nbRang-1),seq(nbRang-1,1))
##aux = length(rang1)
##rang2 <- unlist(lapply(2:nbRang,seq,to=nbRang,by=1))
##rang1 <- rep(rang1,length(jugeActuel))
##rang2 <- rep(rang2,length(jugeActuel))
##jugeActuel <- rep(jugeActuel,each=aux)
##mat <- data.frame(t(matrix(c(rang1,rang2,jugeActuel),ncol=3)))
##
##complete.mat = function(liste, matJR){
##            rang1 <- liste[1]
##            rang2 <- liste[2]
##            jugeActuel <- liste[3]
##            # Sélection de 2 produits à échanger
##            if (matJR[jugeActuel,rang1]<matJR[jugeActuel,rang2])  {
##              prod1<-matJR[jugeActuel,rang1]
##              prod2<-matJR[jugeActuel,rang2]
##            }
##            else {
##              prod1<-matJR[jugeActuel,rang2]
##              prod2<-matJR[jugeActuel,rang1]
##            }
##          return(c(rang1,rang2,jugeActuel,prod1,prod2))  
##            
##}
##mat.comp = lapply(mat, complete.mat,matJR=matJR)
##calc.crit = function(liste, matJR, matPR, unifPR, unifPP, weightRang, weightSucc){
##            rang1 <- liste[1]
##            rang2 <- liste[2]
##            jugeActuel <- liste[3]
##            prod1 <- liste[4]
##            prod2 <- liste[5]
##            # Calcul de la qualité du quadruplet associé dans PxR
##            critPR<-critQualPR(prod1, rang1, prod2, rang2, matPR, unifPR, matJR, jugeActuel)
##            # Calcul de la qualité du quadruplet associé dans PxP
##            critPP<-critQualPP(rang1, rang2, nbRang, matPP, unifPP, matJR, jugeActuel)
##            # Sauvegarde des meilleurs couples...
##            critGlob<-weightRang*critPR+weightSucc*critPP
##return(c(rang1,rang2,jugeActuel,prod1,prod2,critGlob))
##}
##res.crit = lapply(mat.comp, calc.crit,matJR=matJR, matPR=matPR, unifPR=unifPR, unifPP=unifPP, weightRang=weightRang, weightSucc=weightSucc)
##res.crit=as.matrix(as.data.frame(res.crit))
##oo=rev(order(res.crit[6,]))[1]
##critGlob = res.crit[6,oo]
##            if (critGlob>critMax) {
##              critMax<-critGlob
##              iStar<-res.crit[3,oo]
##              prod1Max<-res.crit[4,oo]
##              prod2Max<-res.crit[5,oo]
##              rang1Max<-res.crit[1,oo]
##              rang2Max<-res.crit[2,oo]
##            }
##### Fin modif

      # MAJ de la matrice PxP
      matPP<-majPP(matPP, iStar, matJR, nbRang, rang1Max, rang2Max)
      # MAJ de la matrice PxR
      matPR<-majPR(matPR, rang1Max, rang2Max, iStar, matJR, prod1Max, prod2Max)
      # MAJ de la matrice JxR
      matJR<-majJR(matJR, rang1Max, rang2Max, iStar)

      it<-it+1
    }

    # Calcul du nombre d'écart à l'uniformité PxR
    eamPR<-eam(matPR, unifPR)
    # Calcul de l'EAM associé à PxP
    eamPP<-eam (matPP, unifPP)
    # Calcul de la qualité du plan actuel
    eamPond<-weightRang*eamPR+weightSucc*eamPP

    # Sauvegarde du meilleur plan JxR parmi nbPlan
    if (eamPond < eamPondMin) {
      matJROpt<-matJR
      eamPondMin<-eamPond
    }
  }
  
  return (matJROpt)
}
  
############# Fonctions pour l'algorithme de Fedorov
# Construit un plan optimal (facteur1 x facteur2) en utilisant l'algorithme
# d'échange de Fedorov. Le plan optimal est choisi parmi nbPlan construit.
# @nbMod1     : Nombre de modalités du facteur 1 (juge).
# @nbMod2     : Nombre de modalités du facteur 2 (produit).
# @nbEssPlan  : Nombre d'essais du plan.
# @matCond    : Matrice condensée du plan complet.
# @matEffet   : Matrice des effets du plan complet.
# @tabEssai   : vecteur de numéros d'essais (essais imposés + tirage aléatoire).
# @nbEssImp   : Nombre d'essais imposés.
# @nbPlan     : Nombre de plan à construire.
# @nbIter     : Nombre d'itérations de la procédure d'échange.
# return      : Un plan optimal (facteur1 x facteur2).
# ========================================================================== #
ginv <- function (X, tol = sqrt(.Machine$double.eps)) 
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

fedorov<-function(nbMod1, nbMod2, nbEssPlan, matCond, matEffet, tabEssais, nbEssImp=0, nbDesign=10)
{
  # fonction interne
  # Vérifie qu'un plan est valide en testant si le nombre de produits attribué
  # par juge est correcte.
  # @plan  : Plan à tester.
  # @npj   : Nombre de produits à attribuer à chaque juge.
  # return : Vrai si le plan est valide, faux sinon.

  planValide<-function(plan, npj) {
    test<-unique(rowSums(plan))
    if(length(test)!=1 | test[1]!=npj) return(FALSE)
    else return(TRUE)
  }

    # nbEssaiImp : Nombre d'essais imposés
  # nbEssaiPlan : Nombre d'essais que l'on souhaite conserver
  # nbEssai : Nombre d'essais initiale
  genJPAlea<-function(nbEssaiImp, indEssImp, nbEssaiPlan, nbEssaiPlanComp) {
    res<-vector(mode="integer", length=nbEssaiPlanComp)
    dejaTire<-vector(length=nbEssaiPlanComp)

    if (nbEssaiImp>0) {
      # Recopie des anciennes données
      res[1:nbEssaiImp]<-indEssImp[1:nbEssaiImp]
      # On sauvegarde les numéros d'essais imposés d'indice 1 à nbEssaiImp
      dejaTire[res[1:nbEssaiImp]]<-TRUE
    }

    # Tirage aléatoire ...
    nbTirage<-nbEssaiPlan-nbEssaiImp
    i<-0
    while (i < nbTirage) {
      tirage<-as.integer(1+runif(1)*nbEssaiPlanComp)
      if (dejaTire[tirage]) next

      res[nbEssaiImp+1+i]<-tirage
      dejaTire[tirage]<-TRUE
      i<-i+1
    }
    return (res)
  }

# ========================================================================== #

  tabEssaisCpy <- tabEssais  # Copie de tabEssai
  nbEssTotal<-nbMod1*nbMod2  # Nombre d'essais total
  modTotalEff<-(nbMod1+nbMod2-1)  # Cumul des modalités "avec effet"  =nbMod1-1+nbMod2-1+1
  planOpt<-matrix()  # construction de l'objet résultat
  detOpt <- -9E10   # Déterminant optimal
  iPlan<-0      # compteur de "Plan construit"
  nbMod1ParMod2<-nbEssPlan/nbMod1   # Facteur1 par facteur2
  
    # Cas ou le plan ne peut pas être optimisé
  if (nbEssImp+1 >= nbEssPlan) {
    # Nombre de présence des essais dans le plan
    nbPres<-vector(mode="integer", length=nbEssTotal)
    nbPres[tabEssaisCpy[1:nbEssPlan]]<-nbPres[tabEssaisCpy[1:nbEssPlan]]+1  #! "...<-1" doit être suffisant

    # Construction du résultat
    res<-matrix(0, nrow=nbMod1, ncol=nbMod2)
## Modif Husson
    for (p in 1:nbEssTotal) res[matCond[p,1], matCond[p,2]]<- nbPres[p]
##    for (p in 1:nbEssTotal) {
##      if(nbPres[p]==0) res[matCond[p,1], matCond[p,2]]<- 0
##      else res[matCond[p,1], matCond[p,2]]<-nbPres[p]
##    }
    return(res)
  }


  # Construction de 'nbPlan' plans optimaux puis sélection du meilleur
  while (iPlan < nbDesign) {
    # Extraction des effets des 'nbEssPlan premières essais de 'tabEssai'
    matEffetSub<-matEffet[tabEssaisCpy[1:nbEssPlan], ]
    PX<-t(matEffetSub)%*%matEffetSub  # Calcul de (X'X)
    detPlan<-det(PX, logarithme=FALSE)  # Déterminant de (X'X)

    # Test si la matrice (X'X) est inversible
    if (detPlan<1E-5) {
      tabEssaisCpy<-genJPAlea (nbEssaiImp=nbEssImp, indEssImp=tabEssaisCpy, nbEssaiPlan=nbEssPlan, nbEssaiPlanComp=nbEssTotal) 
      next    # renvoi en debut de boucle while sans incrémentation de iplan
    }
    
    DISP<-ginv(PX)   # Inverse de (X'X)

    # Nombre de présence des essais dans le plan
    nbPres<-vector(mode="integer", length=nbEssTotal)
    nbPres[tabEssaisCpy[1:nbEssPlan]]<-nbPres[tabEssaisCpy[1:nbEssPlan]]+1  #! "...<-1" doit être suffisant
    
    # Procédure d'échange itérative
    while(TRUE) {
      jmax<-0
      kmax<-0
      deltaMax<-0    # Plus grande variation du plan
      # Variance des 'nbEssTotal' essais
      tabVar<-diag((matEffet%*%DISP)%*%t(matEffet))

      # On sort de la boucle si le plan ne peut plus être optimisé
      if (nbEssImp+1 >= nbEssPlan) break
      # sinon on procède aux échanges entre le plan final et le plan complet
      for (k in c((nbEssImp+1):nbEssPlan)) {
        essaiK<-tabEssaisCpy[k]
        covar1<-(matEffet %*% t(DISP))[essaiK, ]  # Calcul premier terme de covariance

        for (j in 1:nbEssTotal) {
          if (nbPres[j]<0.5 && tabVar[j]>tabVar[essaiK]) {
            covkj<-sum(covar1 %*% (t(matEffet))[ ,j])    # Calcul second terme de covariance            
            delta<-tabVar[j]*(1-tabVar[essaiK])-tabVar[essaiK]+covkj^2
            if (delta>deltaMax) {
              deltaMax<-delta
              jmax<-j
              kmax<-k
            }
          }
        }
      }
      
      # si accroissement nul : on passe à l'itération suivante
      if (deltaMax < 1E-15) break

      # sinon MAJ de (X'X)-1 : calcul de la nvelle variance de kmax : kmax.ZZ.kmax
      Z<-(DISP %*% t(matEffet))[ ,jmax]
      ZZ<-matrix(0, modTotalEff, modTotalEff)
      for(i2 in 1:modTotalEff) {
        for(i1 in 1:modTotalEff) ZZ[i1,i2]<-DISP[i1,i2]-((Z[i1]*Z[i2])/(1+tabVar[jmax]))
      }

      matEffetSub<-matEffet[tabEssaisCpy[kmax],]
      vart<-sum((matEffetSub %*% t(ZZ)) %*% matEffetSub)
      Z<-t(ZZ %*% t(matEffet)[ ,tabEssaisCpy[kmax]])

      for (i2 in 1:modTotalEff) {
        for (i1 in 1:modTotalEff) DISP[i1,i2]<-ZZ[i1,i2]+(Z[i1]*Z[i2])/(1-vart)
      }

      # MAJ de nbPres avec les nouveaux paramètres optimaux
      nbPres[jmax]<-nbPres[jmax]+1
      nbPres[tabEssaisCpy[kmax]]<-nbPres[tabEssaisCpy[kmax]]-1
      tabEssaisCpy[kmax]<-jmax

      # Déterminant du plan nouvellement construit
      detPlan<-detPlan*(1+deltaMax)
    } #fin de la procédure d'échange
    

    # Construction du plan juge * produit optimal
    res<-matrix(data=0, nrow=nbMod1, ncol=nbMod2)
## modif Husson
    for(p in 1:nbEssTotal) res[matCond[p,1], matCond[p,2]]<-nbPres[p]
##    for(p in 1:nbEssTotal) {
##      if(nbPres[p]==0) res[matCond[p,1], matCond[p,2]]<-0
##      else res[matCond[p,1], matCond[p,2]]<-nbPres[p]
##    }

    # Test si le plan est valide (utilisé seulement pour les petites matrices
    # car les valeurs d'accroissement des déterminants sont trop faibles pour
    # être détectés correctement.
    #if (nbMod1 <= 5 && nbMod2 <= 5)
      if (!planValide (res, nbMod1ParMod2)) {
        iPlan<-iPlan+1
        if(iPlan >= nbDesign){print(res); stop("\nPas de solution valide trouvée : \nessayez de relancer l'algorithme, ou d'augmenter le nombre de juges ou \nle nombre de produits par juge.\n")}
        else next
      }

    # Sauvegarde du plan ayant le meilleur déterminant
    if (detPlan > detOpt) {
      detOpt <- detPlan
      planOpt <- res
    }

    # Préparation pour la construction du plan suivant...
    tabEssaisCpy<-genJPAlea(nbEssaiImp=nbEssImp, indEssImp=tabEssais, nbEssaiPlan=nbEssPlan, nbEssaiPlanComp=nbEssTotal)  
    iPlan<-iPlan+1
  }
  
  return (planOpt)
}


# ========================================================================== #
# fin des fonctions internes
# ========================================================================== #
  
  # Vérification des paramètres
  if (nbPanelist<nbPanelistMin || nbProdByPanelist>nbProd) stop("CPDO : Paramètres incorrecte !")
  # Initialisation du générateur aléatoire
  set.seed(graine)
  # création de la liste des résultats ($plan, $ordre) à retourner
  res<-list(design=list(), rank=list())
  # Initialisation du Nombre de juges imposé pour l'ordre et de la matrice associée
  nbEssJuge<-0
  matEssJuge<-NA
  nbEssImp <- 0
  if (!is.na(matEssImp)) nbEssImp <- nrow(matEssImp)

  # Construction descendante des plan optimaux de nbPanelistMin à nbPanelist
  for (nbPanelistActuel in (nbPanelistMin : nbPanelist)) {
    nep<-nbPanelistActuel*nbProdByPanelist   # Nombre d'essais du plan
    nep.complet<-nbPanelistActuel*nbProd   # Nombre d'essais du plan complet
    
    # Vérification de l'existance d'une solution
    if (nbPanelistActuel + nbProd - 1 > nep) stop("Matrice non inversible : nombre d'essais insuffisant")

    # Matrice condensée et des effets du plan complet
    matCond<-condenseProd (nbPanelist=nbPanelistActuel, nbProd=nbProd, nbEssai=nep.complet)
    matCandidate<-matEffet(nbPanelist=nbPanelistActuel, nbProd=nbProd, matCond=matCond)

    # Lecture des essais imposés et tirage aléatoire du plan servant à initier Fedorov
    if (nbEssImp > 0) tabEssais<-findEssImp(matEssImp, matCond)
    else tabEssais<-NULL

    if (nbEssImp<nep) tabEssais<-genJPAlea(nbEssaiImp=nbEssImp, indEssImp=tabEssais, nbEssaiPlan=nep, nbEssaiPlanComp=nep.complet)

    # Construction du plan optimal par fedorov
    planOptProd<-fedorov(nbPanelistActuel, nbProd, nep, matCond, matCandidate, tabEssais, nbEssImp, nbDesign = nbDesignProd )
    rownames(planOptProd) = paste("Panelist",1:nrow(planOptProd))
    colnames(planOptProd) = paste("Product",1:ncol(planOptProd))

    # Calcul des essais imposés pour le plan suivant
    lstEssais<-calcEssaisImp (planOptProd)
    nbEssImp<-lstEssais$nb
    matEssImp<-lstEssais$valeurs

    # Ordonnancement des produits
    planOptOrdre<-NA
    if(ordre) {
      planOptOrdre<-rang(nbPanelistActuel, nbProd, nbProdByPanelist, planOptProd,weight, nbPanelistImp = nbEssJuge, matJugeImp = matEssJuge, nbDesignOrdre)
      nbEssJuge<-nrow(planOptOrdre)
      matEssJuge<-planOptOrdre
      rownames(planOptOrdre) = paste("Panelist",1:nrow(planOptOrdre))
      colnames(planOptOrdre) = paste("Rank",1:ncol(planOptOrdre))
    }
    
    # Sauvegarde du plan optimal
##    res$design<-c(res$plan, list(planOptProd))
##    res$rank<-c(res$ordre, list(planOptOrdre))
  }
  res$design<- planOptProd
  res$rank<- planOptOrdre
  class(res)<-list("planDO", "list")
  return (res)
}
