"plotellipseinter" <- function(mat,alpha=0.05,coord=c(1,2),nbgroup=1,moy=TRUE,eig,cex=1,color=NULL,title=NULL,graph.type=c("ggplot","classic")){

#################################################################
"ellipse2" <- function(loc, cov,alpha)
      {
            A <- cov
            detA <- A[1, 1] * A[2, 2] - A[1, 2]^2
            dist <- sqrt(qchisq(1-alpha, 2))
            ylimit <- sqrt(A[2, 2]) * dist
            y <- seq( - ylimit, ylimit, 0.01 * ylimit)
            sqrt.discr <- sqrt(detA/A[2, 2]^2 * abs(A[2, 2] * dist^2 - y^2))
            sqrt.discr[c(1, length(sqrt.discr))] <- 0
            b <- loc[1] + A[1, 2]/A[2, 2] * y
            x1 <- b - sqrt.discr
            x2 <- b + sqrt.discr
            y <- loc[2] + y
            return(rbind(cbind(x1, y), cbind(rev(x2), rev(y))))
      }
#################################################################

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
graph.type <- match.arg(graph.type[1],c("ggplot","classic"))
if (moy ==TRUE) {
  matP=cbind.data.frame(mat$moy$P[,coord],mat$moy$P[,ncol(mat$moy$P)])
  matPJ=cbind.data.frame(mat$moy$PJ[,coord],mat$moy$PJ[,ncol(mat$moy$PJ)])
  matsimul=cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
}
if (moy == FALSE) {
  matmoyP=cbind.data.frame(mat$moy$P[,coord],mat$moy$P[,ncol(mat$moy$P)])
  matmoyPJ=cbind.data.frame(mat$moy$PJ[,coord],mat$moy$PJ[,ncol(mat$moy$PJ)])
  matmoysimul=cbind.data.frame(mat$moy$simul[,coord],mat$moy$simul[,ncol(mat$moy$simul)])
  matP=cbind.data.frame(mat$partiel$P[,coord],mat$partiel$P[,ncol(mat$partiel$P)])
  matPJ=cbind.data.frame(mat$partiel$PJ[,coord],mat$partiel$PJ[,ncol(mat$partiel$PJ)])
  matsimul=cbind.data.frame(mat$partiel$simul[,coord],mat$partiel$simul[,ncol(mat$partiel$simul)])
}  
nbp <- nrow(matP)
nbprod <- nbp/nbgroup
coord.ellipse.a.tracer <- matrix(0,402,2*nbp)

p <- 2
nbjuge <-  nrow(matPJ)/nrow(matP)
nbsimul <-  nrow(matsimul)/nrow(matP)

for (i in 1:nbp){
  VX <- var(matsimul[((i-1)*nbsimul+1):(i*nbsimul),1:2])
  coord.ellipse.a.tracer[,(1+2*(i-1)):(2*i)] <- ellipse2(t(matP[i,1:2]),VX,alpha)
}

minx <- min(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
maxx <- max(coord.ellipse.a.tracer[,1+2*(0:(nbp-1))],na.rm=TRUE)
miny <- min(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)
maxy <- max(coord.ellipse.a.tracer[,2*(1:nbp)],na.rm=TRUE)

lab.x <- paste0("Dim ",coord[1]," (",format(eig[coord[1],2],nsmall=2,digits=2),"%)")
lab.y <- paste0("Dim ",coord[2]," (",format(eig[coord[2],2],nsmall=2,digits=2),"%)")
if (graph.type=="classic"){
  plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = c(minx*1.05,maxx*1.05), ylim = c(1.05*miny,1.05*maxy), col = "white", asp=1)
  if (is.null(title)){
    if (moy==TRUE) title(main = "Confidence ellipses for the mean points")
    if (moy==FALSE) title(main = "Confidence ellipses for the partial points")
  } else {title(main=title)}
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  if (moy==FALSE){
    points(matmoyP[,1], matmoyP[,2],cex=0.8*cex,col=color[1:nbprod],pch=15)
    text( matmoyP[,1], matmoyP[,2], matmoyP[,ncol(matmoyP)], cex = 0.8*cex, pos = 4, offset = 0.2,col=color[1:nbprod])
  }
  if (moy==TRUE) text( matP[,1], matP[,2], matP[,ncol(matP)], cex = 0.8*cex, pos = 4, offset = 0.2,col=color[1:nbprod])
  for (j in 1:nbgroup){
    for (i in 1:nbprod) {
      points(matP[(j-1)*nbprod+i,1], matP[(j-1)*nbprod+i,2],cex=0.8*cex,col=color[i],pch=20)
      if (moy==FALSE) lines(c(matP[(j-1)*nbprod+i,1],matmoyP[i,1]),c( matP[(j-1)*nbprod+i,2],matmoyP[i,2]),col=color[i],lty=j)
      lines(coord.ellipse.a.tracer[,(1+2*((i+(j-1)*nbprod)-1)):(2*(i+(j-1)*nbprod))],col=color[i],lty=j)
    }
  }
} else{
  gg_graph <- ggplot()+ coord_fixed(ratio = 1) +  xlim(c(minx,maxx)) + ylim(c(miny,maxy)) +
    geom_hline(yintercept = 0,lty=2) + geom_vline(xintercept = 0,lty=2) +
    theme_light()+ theme(axis.title = element_text(hjust = 1, face = 2), plot.title = element_text(hjust = 0.5, face = 2))
    if (is.null(title) & moy==TRUE) gg_graph <- gg_graph + labs(x = lab.x, y= lab.y,title = "Confidence ellipses for the mean points")
    if (is.null(title) & moy==FALSE) gg_graph <- gg_graph + labs(x = lab.x, y= lab.y,title = "Confidence ellipses for the partial points")
    if (!is.null(title)) gg_graph <- gg_graph + labs(x = lab.x, y= lab.y,title = title)
    if (moy==FALSE){
      gg_graph <- gg_graph + geom_point(aes(x=matmoyP[,1], y=matmoyP[,2]),color = color[1:nbprod],pch=15,size=0.8*cex*4)
      gg_graph <- gg_graph + ggrepel::geom_text_repel(aes(x=matmoyP[,1], y=matmoyP[,2], label=matmoyP[,ncol(matmoyP)]), color = color[1:nbprod])
    }
    if (moy==TRUE)  gg_graph <- gg_graph + ggrepel::geom_text_repel(aes(x=matP[,1], y=matP[,2], label=matP[,ncol(matP)]), color = color[1:nbprod])
    gg_graph <- gg_graph + geom_point(aes(x=matP[,1], y=matP[,2]),color=rep(color[1:nbprod],nbgroup),pch=20,size=0.8*cex*3)
    for (j in 1:nbgroup){
      for (i in 1:nbprod)  gg_graph <- gg_graph + geom_path(aes_string(x=coord.ellipse.a.tracer[,2*(i+(j-1)*nbprod)-1],y=coord.ellipse.a.tracer[,2*(i+(j-1)*nbprod)]), color =color[i],lty=j)
    }
    if (moy==FALSE){
      group <- as.factor(rep(1:nbgroup,each=nbprod))
      gg_graph <- gg_graph + geom_segment(aes(x=matP[,1],xend=rep(matmoyP[,1],nbgroup),y=matP[,2],yend=rep(matmoyP[,2],nbgroup),lty=group),color=rep(color[1:nbprod],nbgroup)) + theme(legend.position="bottom")
    }
    return(gg_graph)
  }
}
