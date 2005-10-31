"plot.partial" <-
function(coord.partial, H , clab.part = 0, clab.first.part = 0 ,coord=c(1,2)){
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
    fac1 <- factor(rep(1:nbind,nbgroup))
    senso.class(coord.partial[[nbnivo]][[1]],fac1, cpoint = clab.first.part , clabel=0.8, col.moy = rep(nbnivo,length(levels(fac1))),col=rep(nbnivo,length(levels(fac1))),add.plot=TRUE)
    title(sub = "Superimposed representation of the partial clouds" , cex.sub = 1.2, col.sub = "steelblue4", font.sub=3, adj = 0.5, line = 3)
    title(sub = paste("(","comp",coord[1],"-","comp",coord[2],")"), cex.sub = 1.2, col.sub = "steelblue4", font.sub=3, adj = 0.5, line = 4)

}
