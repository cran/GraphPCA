Visu <-
function(PC,Row.names=NULL, labs=NULL,axes=c(1,2),size = 4 )
{
x1=x2=y1=y2=Group=NULL
n=nrow(PC);

x1.var=2*axes[1]-1
x2.var=2*axes[1]
y1.var=2*axes[2]-1
y2.var=2*axes[2]

####################

#auxil<-vector('list',p)
PC1 = PC[,x1.var]
  for(j in c(x2.var,y1.var,y2.var))
      {PC1<-cbind(PC1,PC[,j])
     PC1 = as.data.frame(PC1)
    PC1

}


colnames(PC1) = c('x1','x2','y1','y2')


  if(is.null(Row.names))
{
Row.names=paste(1:n, sep='.')
}
PC1$Group = as.factor(rownames(PC))

if (is.null(labs))
{
labs=paste('Axe',axes,sep='.')
}

plotg = ggplot() + 
scale_x_continuous(name="x") + 
scale_y_continuous(name="y") +
geom_rect(data=PC1, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Group), color="black", alpha=0.5) +labs(x = labs[1], y = labs[2])   

plotg=plotg+geom_text(data=PC1, aes(x2, y2),label = rownames(PC),size=size)


print(plotg)


}
