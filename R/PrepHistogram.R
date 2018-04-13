options(warn=0)
PrepHistogram <-
function(X,Z=NULL,k=3,group=NULL){
hist=NULL

# la fonction ci dessus construit des histogrammes a partir d une
# variable quantitative X et une variable qualitative Z
# Example, si j ai une variable quantitative X=Age et une variable
# qualitative Z=pays, en appelant prephistogram, on construit des
# des histogram pour chaque pays, on n obtient pas des dessins
# mais un tableau dans lequel on a des chiffres pour chaque pays
# et la somme des valeurs de chaque ligne est egal a 1
# le parametre k represente le nombre de classe par histogramme
# les sorties de cette function sont : mido (le centre des classes)
# et Vhistogramme (  le tableau histogramme que nous utiliserons dans 
# lacp des donnees histogrammes



# input
#X: quantitative variable that need to be transform into symbolic histogr (variable quantitative classique dont on veut extraire des histogrammes, age par exemple)
#Z: categorical variable that need to be used for  the clustering the individuals (variable qualitative qui sert de support pour la construction des histogrammes, pays par example)
#k: number of bins of histogram-valued variables (nombre de classes par histogramme, par defaut dans le programme on a k=3 par example)
  

# output
# mido : return the class centres of class and the min and the max of histograms
# histo: dataframe containing the relative frquency of histogram variables
  
n=length(as.data.frame(X))  # on determine n le nombre d observations, pour ce faire on determine la longueur du vecteur colonne X
ncat=nrow(table(Z))  # On determine le nombre de categories de la variable qualitative, si par exemple Z=c("congo","mauritanie"), ncat=2 car il y a deux pays
Vhistogram=matrix(0,ncat,k) # ensuite, on cree une matrice d histogramme de 0, ayant ncat lignes (nombre de categorie de la variables qualitative) et k colonnes (nombre de classe qu on veut assigner a l histogramme)

mido=seq(min(X),max(X),   (range(X)[2]-range(X)[1])/k) # discretisation: on delimite l histogramme

# on va creer une sequence de point allant de min(X) jusqu a max(X), et le pas
# de la sequence est (max(X)-min(X))/2; en fait, range(X)[2]=max(X) et range(X)[1]=min(X)


# on va maintenant remplir Vhistogram ligne par ligne; pour ce faire on va ecrire une 
# soit i allant de 1 a ncat, dans le cas Z=c("congo","mauritanie"), ncat=2 donc i in {1,2}
# soit i=1


for (i in 1:ncat)
{
Vhistogram[i,]<-(hist(subset(X,Z==names(table(Z)[i])),breaks=mido,plot=F)$counts)/(sum(hist(subset(X,Z==names(table(Z)[i])),breaks=mido,plot=F)$counts))
}

## rappel X et Z viennent a la base d une meme base, X est quantitative, Z est qualitative
## subset(X,Z==names(table(Z)[i]) : pour i=1 (index de la premiere categorie de Z c est a dire congo)
## si k=3 par exemple,on va calcule le nombres d observation comprise entre (minX, ....,maxX) , (voir mido l intervalle de discretisation)

# par consequent, Vhistogram[i,]<-(hist(subset(X,Z==names(table(Z)[i])),breaks=mido,plot=F)$counts)/(sum(hist(subset(X,Z==names(table(Z)[i])),breaks=mido,plot=F)$counts)), on fait la meme chose que precedemment, mais 
# on divise par le total

# une fois que c est fait, on reprend la boucle pour i=2

Vhistogram=as.data.frame(Vhistogram)
# on stocke le resultat sous la forme de dataframe

rownames(Vhistogram)=names(table(Z))
# pour chaque categorie de la variable qualitative, on stocke le nom ("congo" et "mauritane")


return(list(mid=mido,Vhistogram=Vhistogram))
# on renvoie "mido" le centre des classes (tres important pour connaitre quelle discretisation est utilisee) et Vhistogram
}
