fast=function(don,alpha=0.05,sep.words=";",word.min=5,graph=TRUE,axes=c(1,2),ncp=5,B=200,val=FALSE, B.val=200,label.miss=NULL){
#don : fichier de données
#alpha : intervalle de confiance pour les ellipses
#sep.words: séparateur de mots
#word.min : fréquence des mots qu'on enlève
#ncp : nombre de composantes principales
#B : nombre de simulations
#val: éléments de validité
#label.miss: label des données incomplètes s'il y en a

don=as.data.frame(don)
I=nrow(don)
J=ncol(don)

#Mise en facteur des variables si ça n'a pas été fait avant
for (i in 1:J){
don[,i]=as.factor(don[,i])}

#ACM
acm=MCA(don,graph=F,ncp=ncp)

if (graph){
plot.MCA(acm,choix="ind",invisible="var",axes=axes) 
plot.MCA(acm,choix="ind",invisible="ind",axes=axes) 
plot.MCA(acm,choix="ind",axes=axes) 
plot.MCA(acm,choix="var",axes=axes)}

########################ENlisib pour les modalités (affiche une partie des modalités)
X=0.05 #Pourcentage de modalités au milieu
tab=as.matrix(acm$var$v.test)                       				  #récupération des vtests des modalités dans un tableau
vtestmean1=mean(tab[,axes[1]])                 					  #calcul de la moyenne des vtests
vtestmean2=mean(tab[,axes[2]])                 					  #calcul de la moyenne des vtests
vtestsd1=sd(tab[,axes[1]])                     					  #calcul de l'écart-type des vtests
vtestsd2=sd(tab[,axes[2]])                     					  #calcul de l'écart-type des vtests
seuil1=vtestmean1+vtestsd1                            				  #seuil=moyenne+écart-type
seuil2=vtestmean2+vtestsd2                            				  #seuil=moyenne+écart-type
modext=which(abs(tab[,axes[1]])>=seuil1 | abs(tab[,axes[2]])>=seuil2)          #récupération des modalités dont la vtest est supérieure au seuil
tab2=which(abs(tab[,axes[1]])<seuil1 & abs(tab[,axes[2]])<seuil2)               #récupération des modalités dont la vtest est inférieure au seuil
modmoy=tab2[sample(1:length(tab2),X*length(tab2))]       					  #tirage au hasard des modalités récupérées à l'étape précédente
mod_kept=cbind(t(modext),t(modmoy))
res=acm                                     					  #remplissage du résultat avec seulement les individus et modalités sélectionnés
res$var$coord=acm$var$coord[mod_kept,]
res$var$cos2=acm$var$cos2[mod_kept,]
res$var$contrib=acm$var$contrib[mod_kept,]
res$var$v.test=acm$var$v.test[mod_kept,]
if (graph){
plot.MCA(res,choix="ind",invisible="ind",axes=axes) }                    #plot individus
########################Fin ENlisib pour les modalités              	

########################Graphiques préliminaires
#Nombre de groupes formés par chaque juge
lev=rep(NA,J)
for (i in 1:J){
lev[i]=nlevels(don[,i])}
lev2=as.factor(lev)
if (graph){
x11()
plot(lev2,main="Number of groups formed from sorting tasks") }

#Nombre de produits par groupe
nb_prod_grp=rep(NA,sum(lev))
grp=0
for (i in 1:J){
nb_prod_grp[(grp+1):(grp+lev[i])]=table(don[,i])
grp=grp+lev[i]}
nb_prod_grp2=as.factor(nb_prod_grp)
if (graph){
x11()
plot(nb_prod_grp2,main="Number of products per group")}
########################Fin graphiques préliminaires


########################Matrice de cooccurrence réordonnée
#Matrice de cooccurrences
tdc=tab.disjonctif(don)
compte=tdc%*%t(tdc)

#Ordre des produits
ordre_prod=order(acm$ind$coord[,1])

#réordonnement matrice cooccurrence
compte2=compte[ordre_prod,]
compte2=compte2[,ordre_prod]
########################Fin matrice de cooccurrence réordonnée

########################Matrice ordonnée des coefficients de Cramer
#Fonction du chi2
chi2_t=function(x,y){
obs=table(x,y)
chi=matrix(NA,length(levels(x)),length(levels(y)))
for (i in 1:length(levels(x))){
for (j in 1:length(levels(y))){
chi[i,j]=(obs[i,j]-(sum(obs[i,])*sum(obs[,j])/length(x)))^2/(sum(obs[i,])*sum(obs[,j])/length(x))
}}
chi2=sum(rowSums(chi))
return(chi2)}

#Fonction du coefficient de Cramer
cramer=function(x, y) {
chi2=chi2_t(x,y)
n=length(x)
p=length(levels(x))
q=length(levels(y))
m=min(p - 1, q - 1)
V=sqrt(chi2/(n * m))
return(V)
}

res=matrix(NA,J,J)
for (i in 1:J){
for (j in i:J){
res[i,j]=res[j,i]=cramer(don[,i],don[,j])}}
colnames(res)=rownames(res)=colnames(don)

afc=CA(res,graph=F)

ord=order(afc$row$coord[,1])
res2=res[ord,]
res2=res2[,ord]
########################Matrice ordonnée des coefficients de Cramer

########################Tableau des données réordonné
out=matrix(NA,I,J)
tdc=tab.disjonctif(don)
gp=0
for (i in 1:J){
conc=cbind(c(1:lev[i]),acm$var$coord[(1+gp):(gp+lev[i]),1])
o=order(conc[,2])
conc2=cbind(conc[o,],c(1:lev[i]))
o2=order(conc2[,1])
conc3=conc2[o2,]
out[,i]=tdc[,(1+gp):(gp+lev[i])]%*%conc3[,3]
gp=gp+lev[i]}

out=data.frame(out)

for (i in 1:J){
out[,i]=as.factor(out[,i])}

#ordre pour les produits
catego_num2=out[ordre_prod,]

#ordre pour les juges
catego_num2=catego_num2[,ord]
rownames(catego_num2)=rownames(don)[ordre_prod]
colnames(catego_num2)=colnames(don)[ord]
########################Tableau des données réordonné


#########################Ellipses
out_axe=function(don,acm){
I=nrow(don)
J=ncol(don)
tdc=tab.disjonctif(don)

axe=vector(mode="list")
axe$eig=acm$eig

vect_prod=rep(rownames(acm$ind$coord),each=J)
vect_juge=rep(1:J,I)
vect_prod_juge=data.frame(vect_prod,vect_juge)
colnames(vect_prod_juge)=c("Product","Panelist")

coord=matrix(NA,(I*J),ncp)
coord=data.frame(coord)
colnames(coord)=paste("Dim",1:ncp,sep=" ")
lev=rep(NA,J)
for (i in 1:J){
lev[i]=length(levels(don[,i]))}
gp=0
for (i in 1:J){
coord[((i-1)*I+1):(i*I),]=tdc[,(1+gp):(gp+lev[i])]%*%acm$var$coord[(1+gp):(gp+lev[i]),1:ncp]
gp=gp+lev[i]}

#Pour avoir les produits au barycentre des mots
for (i in 1:ncp){
coord[,i]=coord[,i]/sqrt(acm$eig[i,1])}

coord.partial=data.frame(coord,vect_prod_juge)

vect_prod_moy=rownames(acm$ind$coord)
vect_juge_moy=rep(0,I)
vect_prod_juge_moy=data.frame(vect_prod_moy,vect_juge_moy)
colnames(vect_prod_juge_moy)=c("Product","Panelist")
coord.partial_moy=data.frame(acm$ind$coord[,1:ncp],vect_prod_juge_moy)

axe$moyen=rbind(coord.partial_moy,coord.partial)
axe$moyen[,"Panelist"]=as.factor(axe$moyen[,"Panelist"])
 
return(axe)}

res.axe=out_axe(don,acm)

sim=simulation(res.axe,nbsimul=B)
if (graph){
x11()
color =  c("black", "red", "green3", "blue", "cyan", "magenta", 
            "darkgray", "darkgoldenrod", "darkgreen", "violet", 
            "turquoise", "orange", "lightpink", "lavender", "yellow", 
            "lightgreen", "lightgrey", "lightblue", "darkkhaki", 
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange", 
            "darkorchid", "darkred", "darksalmon", "darkseagreen", 
            "darkslateblue", "darkslategray", "darkslategrey", 
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon", 
            "lightyellow", "maroon")
repet=I/length(color)
color=rep(color,ceiling(repet))
plotellipse (sim, alpha = alpha, eig = signif(acm$eig,4),coord=axes,color=color)}

#Calcul du rapport
inter=sum(tapply(sim[[1]][[3]][,1],sim[[1]][[3]][,3],mean)^2)/I
tot=sum(sim[[1]][[3]][,1]^2)/(B*I)
ratio=inter/tot
###########################Fin ellipses


###########################Analyse textuelle
texte=matrix(NA,(I*J),3)
texte=data.frame(texte)
texte[,1]=rep(rownames(don),J)
texte[,2]=rep(colnames(don),each=I)
for (i in 1:J){
texte[((I*(i-1))+1):(I*i),3]=paste(don[,i])}
#On ne prend pas le tiret comme séparateur, ni l'apostrophe
restext=textual(texte,3,1,sep.word=sep.words)

#Suppression des modalité g1, ..., g99 (attention tout est mis en minuscule avce textual)
mod.suppr=paste("g",1:99,sep="")
mod.suppr=intersect(colnames(restext$cont.table),mod.suppr)
if (length(mod.suppr)!=0){
num.mod.suppr=which(colnames(restext$cont.table)%in%mod.suppr)
restext$cont.table=restext$cont.table[,-num.mod.suppr]
num.mod.suppr2=which(rownames(restext$nb.words)%in%mod.suppr)
restext$nb.words=restext$nb.words[-num.mod.suppr2,] }

#Nombre de mots différents
nb_mot_diff=nrow(restext$nb.words)
cat("Number of different words : ",nb_mot_diff,"\n")

#Nombre de mots par classe
mots=rep(NA,sum(lev))
grp=0
for (i in 1:J){
mots[(grp+1):(grp+lev[i])]=levels(don[,i])
grp=grp+lev[i]}
mots_split=strsplit(mots,split=sep.words)
nb_mots=rep(NA,length(mots_split))
for (i in 1:length(mots_split)){
if (mots_split[[i]][1] %in% paste("G",1:99,sep="")){
nb_mots[i]=0}
else {
nb_mots[i]=length(mots_split[[i]])}}
nb_mots2=as.factor(nb_mots)

if (graph){
x11()
plot(nb_mots2,main="Number of words per group")}

#Seuil minimum à mettre en paramètre...
freq_min=which(apply(restext$cont.table,2,sum)<=word.min)
if (length(freq_min)!=0){
restext$cont.table=restext$cont.table[,-freq_min]}

juxt=matrix(NA,sum(restext$cont.table),2)
colnames(juxt)=c("Product","Word")
gp=0
for (i in 1:nrow(restext$cont.table)){
for (j in 1:ncol(restext$cont.table)){
if (restext$cont.table[i,j]!=0){
juxt[(gp+1):(gp+restext$cont.table[i,j]),1]=rownames(restext$cont.table)[i]
juxt[(gp+1):(gp+restext$cont.table[i,j]),2]=colnames(restext$cont.table)[j]
gp=gp+restext$cont.table[i,j]}}}
caract_prod=catdes(data.frame(juxt),1)
###########################Analyse textuelle


##################Elements of validity
if (val==TRUE){
eig_perm=ratio_perm=rep(NA,B.val)
for (i in 1:B.val){
don_perm=don
for (j in 1:J){
if (!is.null(label.miss)){
don_perm[which(don[,j]==label.miss),j]=don[which(don[,j]==label.miss),j]
don_perm[which(don[,j]!=label.miss),j]=don[sample(which(don[,j]!=label.miss)),j]}
else {
don_perm[,j]=don[sample(1:I,I),j]}}
acm_perm=MCA(don_perm,graph=F)
eig_perm[i]=acm_perm$eig[1,1]
res.axe_perm=out_axe(don_perm,acm_perm)
sim_perm=simulation(res.axe_perm,nbsimul=B)
inter_perm=sum(tapply(sim_perm[[1]][[3]][,1],sim_perm[[1]][[3]][,3],mean)^2)/I
tot_perm=sum(sim_perm[[1]][[3]][,1]^2)/(B*I)
ratio_perm[i]=inter_perm/tot_perm
}
eig_p.value=length(which(eig_perm>acm$eig[1,1]))/B.val
ratio_p.value=length(which(ratio_perm>ratio))/B.val
}
##################End elements of validity

##################Sorties
acm_call=list(X=acm$call$X,marge.col=acm$call$marge.col,marge.row=acm$call$marge.row,ncp=acm$call$ncp,quali=acm$call$quali,mca=acm,sim=sim)
group_afm=list(coord=acm$var$eta2)
if(val==TRUE){
validity=vector(mode="list")
validity$eig=vector(mode="list")
validity$ratio=vector(mode="list")
validity$eig$value=acm$eig[1,1]
validity$eig$p.value=eig_p.value
validity$eig$permut=eig_perm
validity$ratio$value=ratio
validity$ratio$p.value=ratio_p.value
validity$ratio$permut=ratio_perm}
else {
validity=NULL}
res=list(eig=acm$eig,var=acm$var,ind=acm$ind,group=group_afm,acm=acm,call=acm_call,cooccur=compte2,reord=catego_num2,cramer=res2,textual=caract_prod,validity=validity)                                                       
class(res) <- c("catego", "list ")
return(res)}
