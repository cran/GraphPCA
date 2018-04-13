Ridi <-
function(X){
. <- "Shut up"
X=as.data.frame(X)
#X=scale(X,scale=F)
entree=apply(X,2,mean)

p=ncol(X)
#res=vector("list",p)
res=entree
res[1]=0.5*entree[1]



for(j in 2:p){
#if (j>1){ res[j]=entree[j]+entree[j-1]}
aux=entree[j-1]
res[j]=0.5*entree[j]+max(cumsum(entree[1:(j-1)]))
}



return(list(p=p,entree=entree,Ridit=res))
}
