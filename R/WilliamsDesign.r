WilliamsDesign = function (nbprod,seed=NULL) {
  if (nbprod == 2) plan <- matrix(c(1, 2, 2, 1), 2, 2)
  if (nbprod > 2) {    
    plan <- matrix(0,nbprod,nbprod)
    ligne = c(0,rep(1,nbprod-1))
    ligne[(2:(nbprod%/%2))*2] = 2:(nbprod%/%2)
    ligne[(2:(nbprod%/%2+nbprod%%2))*2-1] = (nbprod-1):((nbprod%/%2+1))
    plan[1,]=ligne
    for (i in 2:nbprod)  plan[i,] <- (plan[i-1,]+1)%%nbprod
    if ((nbprod%%2)==1) {
      plan2 <- plan
      plan2[1,] <- rev(ligne)
      for (i in 2:nbprod)  plan2[i,] <- (plan2[i-1,]+1)%%nbprod
      plan <- rbind(plan,plan2)
    }
    row.names(plan) <- NULL
    plan <- plan+1
  }
  if (!is.null(seed)) {
    set.seed(seed)
    alea <- nbprod+sample(1:nbprod)
    for (k in 1:nbprod) plan =matrix(replace(plan,plan==k,alea[k]),ncol=nbprod)
    plan <- plan-nbprod
  }
  return(plan)
}
