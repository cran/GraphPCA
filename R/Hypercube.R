Hypercube <-
function(g)
{
g<-as.data.frame(g)
l<-nrow(g)
out<-vector('list',l)

for(i in 1:l)
{
out[[i]]<-cube(g[i,])
}
result<-as.data.frame(out[[1]])
for(i in 2:l)
{
res<-as.data.frame(out[[i]])
result<-rbind(result,res)
}
return(list(result=result))
}
