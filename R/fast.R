fast=function(don,alpha=0.05,mot_min=2,graph=TRUE,ncp=5,B=200){
#don : fichier de données
#alpha : intervalle de confiance pour les ellipses
#mot_min : fréquence des mots qu'on enlève
#ncp : nombre de composantes principales
#B : nombre de simulations

don=as.data.frame(don)
I=nrow(don)
J=ncol(don)

#Mise en facteur des variables si ça n'a pas été fait avant
for (i in 1:J){
don[,i]=as.factor(don[,i])}

#ACM
acm=MCA(don,graph=F,ncp=ncp)
#rownames(acm$var$coord)=colnames(tab.disjonctif(don))
if (graph){
plot.MCA(acm,choix="ind",invisible="var") 
plot.MCA(acm,choix="ind",invisible="ind") 
plot.MCA(acm,choix="ind") }

#ENlisib pour les modalités
X=0.05 #Pourcentage de modalités au milieu
tab=as.matrix(acm$var$v.test)                       				  #récupération des vtests des modalités dans un tableau
vtestmean1=mean(tab[,1])                 					  #calcul de la moyenne des vtests
vtestmean2=mean(tab[,2])                 					  #calcul de la moyenne des vtests
vtestsd1=sd(tab[,1])                     					  #calcul de l'écart-type des vtests
vtestsd2=sd(tab[,2])                     					  #calcul de l'écart-type des vtests
seuil1=vtestmean1+vtestsd1                            				  #seuil=moyenne+écart-type
seuil2=vtestmean2+vtestsd2                            				  #seuil=moyenne+écart-type
modext=which(abs(tab[,1])>=seuil1 | abs(tab[,2])>=seuil2)          #récupération des modalités dont la vtest est supérieure au seuil
tab2=which(abs(tab[,1])<seuil1 & abs(tab[,2])<seuil2)               #récupération des modalités dont la vtest est inférieure au seuil
modmoy=tab2[sample(1:length(tab2),X*length(tab2))]       					  #tirage au hasard des modalités récupérées à l'étape précédente
mod_kept=cbind(t(modext),t(modmoy))
res=acm                                     					  #remplissage du résultat avec seulement les individus et modalités sélectionnés
res$var$coord=acm$var$coord[mod_kept,]
res$var$cos2=acm$var$cos2[mod_kept,]
res$var$contrib=acm$var$contrib[mod_kept,]
res$var$v.test=acm$var$v.test[mod_kept,]
if (graph){
plot.MCA(res,choix="ind",invisible="ind") }                    #plot individus
              	


#Matrice de cooccurrence
compte=matrix(0,I,I)
for (i in 1:J){
for (j in 1:I){
for (k in j:I){
if (don[j,i]==don[k,i]){
compte[j,k]=compte[j,k]+1}}}}

for (i in 1:I){
for (j in i:I){
compte[j,i]=compte[i,j]}}
colnames(compte)=rownames(compte)=rownames(don)


#Matrice des coefficients de Cramer
chi2_t=function(x,y){
obs=table(x,y)
chi=matrix(NA,length(levels(x)),length(levels(y)))
for (i in 1:length(levels(x))){
for (j in 1:length(levels(y))){
chi[i,j]=(obs[i,j]-(sum(obs[i,])*sum(obs[,j])/length(x)))^2/(sum(obs[i,])*sum(obs[,j])/length(x))
}}
chi2=sum(rowSums(chi))
return(chi2)}

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


#Réordonnement de la matrice des données
afm=MFA(don,group=rep(1,J),type=rep("n",J),graph=F,name.group=colnames(don))

if (graph){
plot.MFA(afm,choix="group")}

ordre_prod=order(acm$ind$coord[,1])

#changement des numérotation
out=matrix(NA,I,J)
lev=rep(NA,J)
for (i in 1:J){
lev[i]=length(levels(don[,i]))}
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

#réordonnement matrice cooccurrence
compte2=compte[ordre_prod,]
compte2=compte2[,ordre_prod]


lev2=as.factor(lev)
if (graph){
x11()
plot(lev2,main="Number of groups formed from sorting tasks") }

#Ellipses
out1=matrix(NA,(I*J),4)
lev=rep(NA,J)
for (i in 1:J){
lev[i]=length(levels(don[,i]))}
tdc=tab.disjonctif(don)
gp=0
for (i in 1:J){
out1[((i-1)*I+1):(i*I),1:2]=tdc[,(1+gp):(gp+lev[i])]%*%acm$var$coord[(1+gp):(gp+lev[i]),1:2]
gp=gp+lev[i]}

#Répétition des produits
prod=rep(rownames(don),J)
#Répétition des juges
jug=rep(colnames(don),each=I)

out2=data.frame(out1)
out2[,3]=prod
out2[,4]=jug
out2[,3]=as.factor(out2[,3])
out2[,4]=as.factor(out2[,4])

d12=acm$ind$coord[,1:2]
pr=rownames(don)
ju=rep("0",I)
out22=data.frame(d12,pr,ju)
colnames(out22)=colnames(out2)
out3=rbind(out22,out2)
out3[,3]=as.factor(out3[,3])
out3[,4]=as.factor(out3[,4])
out4=out3
out4[-c(1:I),1]=out3[-c(1:I),1]/sqrt(acm$eig[1,1])
out4[-c(1:I),2]=out3[-c(1:I),2]/sqrt(acm$eig[2,1])
out5=out4
out5$moyen=out5

sim=simulation(out5,nbsimul=B)
if (graph){
x11()
plotellipse (sim, alpha = alpha, eig = signif(acm$eig,4))}


#Analyse textuelle
texte=matrix(NA,(I*J),3)
texte=data.frame(texte)
texte[,1]=rep(rownames(don),J)
texte[,2]=rep(colnames(don),each=I)
for (i in 1:J){
texte[((I*(i-1))+1):(I*i),3]=paste(don[,i])}
#On ne prend pas le tiret comme séparateur, ni l'apostrophe
#restext=textual(texte,3,1,sep.word="(), ?;/:'!$=+\n;{}<>[]-")
#restext=textual(texte,3,1,sep.word="(), ?;/:!$=+\n;{}<>[]")

#Nombre de mots différents
#nb_mot_diff=nrow(restext$nb.words)
#cat("Number of different words : ",nb_mot_diff,"\n")

#Nombre de parfums par classe
nbp=strsplit(summary(don,maxsum=max(lev)),":")
agg=rep(0,J*max(lev))
for (i in 1:(J*max(lev))){
agg[i]=nbp[[i]][2]}
agg2=na.omit(agg)
agg2=as.factor(agg2)
if (graph){
x11()
plot(agg2,main="Number of products per group")}

#Nombre de mots par classe
#nb_mot=rep(0,J*max(lev))
#for (i in 1:(J*max(lev))){
#nb_mot[i]=nbp[[i]][1]}
#nb_mot2=na.omit(nb_mot)
#nb_mot2=as.factor(nb_mot2)
#nb_mot3=as.character(nb_mot2)
#nb_mot4=chartr("123456789", "444444444", nb_mot3)
#nb_m=rep(NA,length(nb_mot4))
#for (i in 1:length(nb_mot4)){
#nb_m[i]=length(which((strsplit(nb_mot4," ")[i][[1]]!="")&((strsplit(nb_mot4," ")[i][[1]]!="Gr4"))))}
#nb_m=as.factor(nb_m)
#
##Seuil minimum à mettre en paramètre...
#freq_min=which(apply(restext$cont.table,2,sum)<=mot_min)
#if (length(freq_min)!=0){
#restext$cont.table=restext$cont.table[,-freq_min]}
#

#juxt=matrix(NA,sum(restext$cont.table),2)
#gp=0
#for (i in 1:nrow(restext$cont.table)){
#for (j in 1:ncol(restext$cont.table)){
#if (restext$cont.table[i,j]!=0){
#juxt[(gp+1):(gp+restext$cont.table[i,j]),1]=rownames(restext$cont.table)[i]
#juxt[(gp+1):(gp+restext$cont.table[i,j]),2]=colnames(restext$cont.table)[j]
#gp=gp+restext$cont.table[i,j]}}}
#restexttot=catdes(data.frame(juxt),1)


#Résultats en sortie
acm_call=list(X=acm$call$X,marge.col=acm$call$marge.col,marge.row=acm$call$marge.row,ncp=acm$call$ncp,quali=acm$call$quali,name.group=afm$call$name.group,sim=sim)
group_afm=list(coord=afm$group$coord,cos2=afm$group$cos2,contrib=afm$group$contrib)
res=list(eig=acm$eig,var=acm$var,ind=acm$ind,group=group_afm,call=acm_call,cooccur=compte2,reord=catego_num2,cramer=res2)                                                       
class(res) <- c("catego", "list ")
return(res)}

