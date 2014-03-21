auxilary <-
function(Tab)
{
p2<-(ncol(Tab));
p<-p2/2;
res<-vector('list',p)

m<-seq(1,p2,2)
n<-seq(2,p2,2)

for(i in 1:p)
{

res[[i]]<-as.numeric( Tab[m[i]:n[i]] );


}
return((res))

}
