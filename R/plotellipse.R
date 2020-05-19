"plotellipse" <- function(mat,alpha=0.05,coord=c(1,2),eig,cex=1,color=NULL,title=NULL,graph.type=c("ggplot","classic")){

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
graph.type <- match.arg(graph.type[1],c("ggplot","classic"))
gg_graph <- plotellipseinter(mat,alpha=alpha,coord=coord,nbgroup=1,moy=TRUE,eig=eig,color=color,cex=cex,title=title,graph.type=graph.type)
res <- list()
res$plotIndEll <- gg_graph
if (graph.type=="ggplot") print(gg_graph)
 if (length(mat$partiel)!=0) {
   if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()
   nbgroup <- length(mat$partiel$simul[,ncol(mat$partiel$simul)]) / length(mat$moy$simul[,ncol(mat$moy$simul)])
   g2 <- plotellipseinter(mat,alpha=alpha,coord=coord,nbgroup=nbgroup,moy=FALSE,eig=eig,color=color,cex=cex,title=title,graph.type=graph.type)
   if (graph.type=="ggplot") print(g2)
   res$plotIndEllPar <- g2
 }
if (graph.type=="ggplot") return(res)
}
