Visu <-
function(PC,Row.names=NULL, labs=NULL,axes=c(1,2) )
{

n=nrow(PC);


x1.var=2*axes[1]-1
x2.var=2*axes[1]
y1.var=2*axes[2]-1
y2.var=2*axes[2]

 x1 <- names(PC)[x1.var]
   x2 <- names(PC)[x2.var]
     y1 <- names(PC)[y1.var]
      y2 <- names(PC)[y2.var]


  if(is.null(Row.names))
{
Row.names=paste(1:n, sep='.')
}

if (is.null(labs))
{
labs=paste('Axe',axes,sep='.')
}

ltext <- Row.names

plotg=ggplot(PC, aes_string(x = x1, y = y1)) + geom_point()

plotg=plotg+geom_rect(data = PC, aes_string(xmin = x1, xmax = x2,  ymin = y1, ymax =y2),fill =  alpha(1:n, 
                2/3))+labs(x = labs[1], y = labs[2])

plotg=plotg+geom_text(aes_string(x2, y2),label = ltext,pch=1:n )

PC$z1=(PC[,x1.var] + PC[,x2.var])/2
PC$z2=(PC[,y1.var] + PC[,y2.var])/2


print(plotg)


}
