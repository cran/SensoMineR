"plotpanelist" <- function(mat,coord=c(1,2),name=FALSE,eig,cex=1,color=NULL,graph.type=c("ggplot","classic"))   {

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
graph.type <- match.arg(graph.type[1],c("ggplot","classic"))
mat <- cbind.data.frame(mat[,coord],mat[,(ncol(mat)-1):ncol(mat)])
nbprod <- length(levels(mat[,ncol(mat)-1]))
lab <- mat[1:nbprod,ncol(mat)-1]
nbjuge <-  length(levels(mat[,ncol(mat)]))-1
nbpoint=nbprod*(nbjuge+1)

minx <- min(mat[1:nbpoint,1],na.rm=TRUE)
maxx <- max(mat[1:nbpoint,1],na.rm=TRUE)
miny <- min(mat[1:nbpoint,2],na.rm=TRUE)
maxy <- max(mat[1:nbpoint,2],na.rm=TRUE)
if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) dev.new()

if (graph.type=="classic"){
 plot(0, 0, xlab = paste("Dim ",coord[1]," (",eig[coord[1],2],"%)",sep=""), ylab = paste("Dim ",coord[2]," (",eig[coord[2],2],"%)",sep=""), xlim = c(minx*1.05,maxx*1.05), ylim = c(1.05*miny,1.05*maxy), col = "white", asp=1)
   abline(v=0,lty=2)
   abline(h=0,lty=2)
   points(mat[1:nbprod,1],mat[1:nbprod,2],cex=cex*1.2,col=color[1:nbprod],pch=15)
   text( mat[1:nbprod,1], mat[1:nbprod,2], labels=lab, cex = cex*0.8, pos = 4, offset = 0.2,col=color[1:nbprod])
   if (name==FALSE) for (i in 1:nbjuge) points(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],cex=cex*0.8,col=color[1:nbprod],pch=20)
   if (name==TRUE) for (i in 1:nbjuge) text(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],labels=mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),ncol(mat)],cex=cex*0.6,col=color[1:nbprod],pch=20)
   title(main = "Individual description")
} else{ 
  gg_graph <- ggplot() +
  coord_fixed(ratio = 1) +
  xlim(c(minx,maxx)) + ylim(c(miny,maxy)) +
  geom_hline(yintercept = 0,lty=2) + geom_vline(xintercept = 0,lty=2) +
  theme_light()+ theme(axis.title = element_text(hjust = 1, face = 2), plot.title = element_text(hjust = 0.5, face = 2))+
  labs(title = "Individual description", x = paste0("Dim ",coord[1]," (",eig[coord[1],2],"%)"), y= paste0("Dim ",coord[2]," (",eig[coord[2],2],"%)"))+
  geom_point(aes(x=mat[1:nbprod,1],y=mat[1:nbprod,2]),col=color[1:nbprod],pch=15)+
  ggrepel::geom_text_repel(aes(x=mat[1:nbprod,1], y=mat[1:nbprod,2], label=rownames(mat)[1:nbprod]), color = color[1:nbprod])+
  geom_point(aes(x=mat[-(1:nbprod),1],y=mat[-(1:nbprod),2]),col=rep(color[1:nbprod],nrow(mat)/nbprod-1),pch=20)+
  if (name==TRUE) geom_text(aes(x=mat[-(1:nbprod),1],y=mat[-(1:nbprod),2],label=mat[-(1:nbprod),4]),col=rep(color[1:nbprod],nrow(mat)/nbprod-1),size=3)
  print(gg_graph)
  return(gg_graph)
  }
# if (name==FALSE) for (i in 1:nbjuge) points(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],cex=cex*0.8,col=color[1:nbprod],pch=20)
# if (name==TRUE) for (i in 1:nbjuge) text(mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),1],mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),2],labels=mat[(nbprod+1+nbprod*(i-1)):(nbprod+nbprod*i),ncol(mat)],cex=cex*0.6,col=color[1:nbprod],pch=20)
}
