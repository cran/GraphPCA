cube <-
function(g)
{
n<-nrow(g)
output<-vector('list',n)

for(i in 1:n)
{
output[[i]]<-expand.grid(auxilary(g[i,]))
}
return(output)
}
