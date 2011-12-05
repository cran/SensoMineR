####################################Function for hierarchical sorting task
fahst=function(don,group,alpha=0.05,graph=TRUE,axes=c(1,2),name.group=NULL,ncp=5,B=100,val=FALSE, B.val=200){
#don: jeu de données
#group: vecteur du nombre de hiérarchie pour chaque sujet
#alpha: intervalle de confiance pour les ellipses
#graph: graphique ou non en sortie
#axes: dimensions représentées sur les graphs
#ncp: nombre de composantes principales
#B: nombre de rééchantillonnages pour les ellipses
#val: éléments de validité
#B.val: nombre de réechantillonnage pour les éléments de validité

don=as.data.frame(don)
I=nrow(don)
J=length(group)

#Mise en facteur des variables si ça n'a pas été fait avant
for (i in 1:ncol(don)){
don[,i]=as.factor(don[,i])}

if (is.null(name.group)) 
name.group <- paste("Sj", 1:length(group), sep = ".")

########################################Preliminary graphs
if (graph){

#Nombre de niveaux par sujet
group2=as.factor(group)
plot(group2,main="Number of levels per subject")

#Nombre de groupes au niveau 1 de chaque sujet
niv1=cumsum(group)-group+1
lev1=rep(NA,J)
for (i in 1:J){
lev1[i]=length(levels(don[,niv1[i]]))}
lev1b=as.factor(lev1)
dev.new()  
plot(lev1b,main="Number of groups formed from first levels")

#Nombre de groupes au dernier niveau de chaque sujet
nivd=cumsum(group)
levd=rep(NA,J)
for (i in 1:J){
levd[i]=length(levels(don[,nivd[i]]))}
levdb=as.factor(levd)
dev.new()
plot(levdb,main="Number of groups formed from last levels")

#Nombre d'objets par groupe au niveau 1 de chaque sujet
nbp1=strsplit(summary(don[,niv1],maxsum=max(lev1)),":")
agg1=rep(0,J*max(lev1))
for (i in 1:(J*max(lev1))){
agg1[i]=nbp1[[i]][2]}
agg1b=na.omit(agg1)
agg1b=as.factor(agg1b)
dev.new()
plot(agg1b,main="Number of objects per group for the first levels")

#Nombre d'objets par groupe au dernier niveau de chaque sujet
nbpd=strsplit(summary(don[,nivd],maxsum=max(levd)),":")
aggd=rep(0,J*max(levd))
for (i in 1:(J*max(levd))){
aggd[i]=nbpd[[i]][2]}
aggdb=na.omit(aggd)
aggdb=as.factor(aggdb)
dev.new()
plot(aggdb,main="Number of objects per group for the last levels")

}
########################################Fin graphiques préliminaires

#AFM
afm=MFA(don,group=group,type=rep("n",J),name.group=name.group,graph=F,ncp=ncp)


##################rapport de corrélation##################################
eta2 <- function(x, gpe) {
vartot <- function(x) {
res <- sum((x - mean(x))^2)
return(res)}
varinter <- function(x, gpe) {
moyennes <- tapply(x, gpe, mean)
effectifs <- tapply(x, gpe, length)
res <- (sum(effectifs * (moyennes - mean(x))^2))
return(res)}
res <- varinter(x, gpe)/vartot(x)
return(res)}
###########################################################################

#calcul pour chaque dimension du rapport de corrélation
coord.niv=matrix(0,ncol(don),ncp)
rownames(coord.niv)=colnames(don)
colnames(coord.niv)=colnames(afm$ind$coord)
for (i in 1:ncol(don)){
coord.niv[i,]=apply(afm$ind$coord,2,eta2,don[,i])}

########################################Graphiques de l'AFM
if (graph){
#Graph des individus
plot.MFA(afm,choix="ind",invisible="quali",axes=axes)
#Graph des mots
plot.MFA(afm,choix="ind",invisible="ind",axes=axes)
#Graph des individus et des mots
plot.MFA(afm,choix="ind",axes=axes)
#Graph des groupes au niveau des sujets
plot.MFA(afm,choix="group",axes=axes)
#Graph des groupes au niveau des sujets et des niveaux
#plot.MFA(afm,choix="group",habillage="group",axes=axes)
#points(coord.niv[,axes],pch=2,col=rep(2:(J+1),times=group))
#text(coord.niv[,axes],label=rownames(coord.niv),col=rep(2:(J+1),times=group),pos=3)
#Graph des niveaux
dev.new(width = 8, height = 8)
plot(coord.niv[,axes], xlab = paste("Dim ", axes[1], " (", signif(afm$eig[axes[1],
2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(afm$eig[axes[2],
2], 4), "%)", sep = ""), xlim = c(0,1), ylim = c(0, 1), pch = 17, main = "Levels representation")
text(coord.niv[, axes[1]], y = coord.niv[, axes[2]], labels = rownames(coord.niv),pos = 3)
#Graph des niveaux et trajectoires
dev.new(width = 8, height = 8)
plot(coord.niv[,axes], xlab = paste("Dim ", axes[1], " (", signif(afm$eig[axes[1],
2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(afm$eig[axes[2],
2], 4), "%)", sep = ""), xlim = c(0,1), ylim = c(0, 1), pch = 17, main = "Levels representation and trajectories")
text(coord.niv[, axes[1]], y = coord.niv[, axes[2]], labels = rownames(coord.niv),pos = 3)
subj=0 
for (j in 1:length(group)){
if (group[j]!=1){
for (i in 1:(group[j]-1)){
lines(x=coord.niv[(subj+i):(subj+i+1),axes[1]],y=coord.niv[(subj+i):(subj+i+1),axes[2]])}}
subj=subj+group[j]
}
}
########################################Fin graphiques de l'AFM

########################################Ellipses
out_axe_mfa=function(don,mfa){
I=nrow(don)
J=nrow(mfa$group$coord)

axe=vector(mode="list")
axe$eig=mfa$eig

vect_prod=rep(rownames(mfa$ind$coord),each=J)
vect_juge=rep(1:J,I)
vect_prod_juge=data.frame(vect_prod,vect_juge)
colnames(vect_prod_juge)=c("Product","Panelist")
coord.partial=data.frame(mfa$ind$coord.partiel[,1:ncp],vect_prod_juge)

vect_prod_moy=rownames(mfa$ind$coord)
vect_juge_moy=rep(0,I)
vect_prod_juge_moy=data.frame(vect_prod_moy,vect_juge_moy)
colnames(vect_prod_juge_moy)=c("Product","Panelist")
coord.partial_moy=data.frame(mfa$ind$coord[,1:ncp],vect_prod_juge_moy)

axe$moyen=rbind(coord.partial_moy,coord.partial)
axe$moyen[,"Panelist"]=as.factor(axe$moyen[,"Panelist"])

axe$moyen=axe$moyen[order(axe$moyen[,ncol(axe$moyen)]),]

return(axe)}

res.axe=out_axe_mfa(don,afm)
simul <- simulation(res.axe,nbsimul =B)
if (graph){
dev.new()
plotellipse (simul, alpha = alpha,coord =axes, eig = signif(res.axe$eig,4))}

#Calcul du rapport sur l'axe 1
inter=sum(tapply(simul$moy$simul[,1],simul$moy$simul[,(ncp+1)],mean)^2)/I
tot=sum(simul$moy$simul[,1]^2)/(B*I)
ratio=inter/tot
########################################Ellipses

##################Elements of validity
if (val==TRUE){
print("The procedure for the elements of validity can be time-consuming")
eig_perm=ratio_perm=rep(NA,B.val)
for (i in 1:B.val){
don_perm=don
num_juge=0
for (j in 1:J){
don_perm[,(num_juge+1):(num_juge+group[j])]=don[sample(1:I,I),(num_juge+1):(num_juge+group[j])]
num_juge=num_juge+group[j]}
mfa_perm=MFA(don_perm,group=group,type=rep("n",J),graph=FALSE)
eig_perm[i]=mfa_perm$eig[1,1]

res.axe_perm=out_axe_mfa(don_perm,mfa_perm)
sim_perm=simulation(res.axe_perm,nbsimul=B)
inter_perm=sum(tapply(sim_perm$moy$simul[,1],sim_perm$moy$simul[,(ncp+1)],mean)^2)/I
tot_perm=sum(sim_perm$moy$simul[,1]^2)/(B*I)
ratio_perm[i]=inter_perm/tot_perm
}
eig_p.value=length(which(eig_perm>afm$eig[1,1]))/B.val
ratio_p.value=length(which(ratio_perm>ratio))/B.val
}
##################End elements of validity

if(val==TRUE){
validity=vector(mode="list")
validity$eig=vector(mode="list")
validity$ratio=vector(mode="list")
validity$eig$value=afm$eig[1,1]
validity$eig$p.value=eig_p.value
validity$eig$permut=eig_perm
validity$ratio$value=ratio
validity$ratio$p.value=ratio_p.value
validity$ratio$permut=ratio_perm}
else {
validity=NULL}

afm_call = list(X = afm$call$X, col.w = afm$call$col.w, 
        row.w = afm$call$row.w, ncp = afm$call$ncp, 
        name.group = afm$call$name.group, simul = simul,group=group,mfa=afm)
var_afm = list(coord=afm$quali.var$coord,contrib=afm$quali.var$contrib,cos2=afm$quali.var$cos2,v.test=afm$quali.var$v.test,coord.lev=coord.niv)

res = list(eig = afm$eig, var = var_afm, ind = afm$ind, group = afm$group, validity=validity,
        call = afm_call)
class(res) <- c("fahst", "list ")

return(res)}