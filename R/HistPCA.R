
HistPCA <-
function(Variable=list,score=NULL,t=1.1, axes=c(1,2),Row.names=NULL,xlim=NULL,ylim=NULL,xlegend=NULL,ylegend=NULL,Col.names=NULL,transformation=1,
                  method='hypercube',proc=0,plot3d.table = NULL,axes2 = c(1, 2, 3), ggplot=1)
{

plot=text=abline=points=rect=cor=dev.new=segments=par=NULL

fintab <- function(gen1) {
        k = nrow(gen1)
        k3 = k/2
        k1 = 2 * ncol(gen1)
        k2 = k1/2
        t5 = vector("list", k2)
        t1 = gen1[1:(k/2), ]
        t1 = data.frame(t1)
        t2 = gen1[((k/2) + 1):k, ]
        t2 = data.frame(t2)
        t3 = cbind(t1, t2)
        t3 = data.frame(t3)
        for (i in 1:k2) t5[[i]] = cbind(t3[, i], t3[, i + k2])
        t6 = vector("list", k2)
        for (i in 1:k2) t6[[i]] = matrix(0, nrow = k3, ncol = 2)
        for (j in 1:k3) {
            for (i in 1:k2) t6[[i]][j, ] = t(as.data.frame(range((t5[[i]])[j, 
                ])))
        }
        t7 = t6[[1]]
        {
            for (s in 2:k2) t7 = cbind(t7, t6[[s]])
        }
        t7 = data.frame(t7)
        return(t7)
    }





f <- function(dimension) {
rnorm=dist=lm.fit=plot=segments=NULL
   rorthogonal <- function(n = 2, dimension) {
       # generate  random orthonormal vectors
      # in dimension dimensional space
      x <- matrix(rnorm(n*dimension), nrow=n, ncol=dimension)
       t(eigen(crossprod(x))$vectors)
  }
  cube <- as.matrix(do.call(expand.grid, rep(list(c(-1,1)),dimension)))
  # the points pairs with distance 2 are connected by edges
  ij <- which(abs(as.matrix(dist(cube))-2)<0.00091188196, arr.ind=TRUE)
  # only need edge in one direction
  ij <- ij[ij[,1]<ij[,2],]
  basis <- rorthogonal(n=2, dim=dimension)
   # p is projection of vectices to subspace of basis
   p <- t(lm.fit(x=t(basis), y=t(cube))$coef)
   # plot vertices
   plot(p[,1],p[,2],asp=1)
  # plot edges
   segments(p[,1][ij[,1]], p[,2][ij[,1]], p[,1][ij[,2]], p[,2][ij[,2]])
return(list(p=p,ij=ij,cube=cube,basis=basis))
 }

#######################################################################


auxilary<-function(Tab)
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


#merci<-auxilary(rec)

cube<-function(g)
{
n<-nrow(g)
output<-vector('list',n)

for(i in 1:n)
{
output[[i]]<-expand.grid(auxilary(g[i,]))
}
return(output)
}



#expand.grid(auxilary(Expert2[1,-1]))


#expand.grid((merci))

#cube(Expert1[,-1])

#expand.grid(list(as.numeric(rec[1:2]),as.numeric(rec[3:4]),as.numeric(rec[5:6])))


Hypercube<-function(g)
{
#construire les hypercube de chaque tableau;
g<-as.data.frame(g)
l<-nrow(g)




out<-vector('list',l)

for(i in 1:l)
{
out[[i]]<-cube(g[i,])
}




#empiler les hypercube par un rbind;
result<-as.data.frame(out[[1]])




for(i in 2:l)
{
res<-as.data.frame(out[[i]])
result<-rbind(result,res)
}
return(list(result=result))
#return(result)
#sortir le tableau pour lancer statis
}


  
  ucad <- function(p) {
    res <- vector("list", p)
    auxi <- NULL
    for (i in 1:p) {
      res[[i]] <- (paste(c("PCMin", "PCMax"), i, sep = "."))
      auxi <- c(auxi, res[[i]])
    }
    auxi <- (as.data.frame(auxi))
    auxi
  }
  
  Centrage <- function(x, reduire = 0) {
    x1 <- scale(x, center = TRUE, scale = FALSE)
    x2 <- scale(x, center = TRUE, scale = TRUE)
    if (reduire == 1) {
      y = x1
    }
    else if (reduire == 0) {
      y = x2
    }
    return(as.matrix(y))
  }
  
  
  Build3d <- function(dimens, GenPC, p, axes2, n, Row.names, axes, 
                      lab.x, lab.y, nametable) {
					  #bibi
					  dev.new=text=plot=abline=rect=NULL
    dev.new()
    k = dimens
    PCinterval1 <- fintab(as.data.frame(GenPC[[k]]))
    tempGenPC <- GenPC[[k]]
    regad <- ucad(p)
    p2 = nrow(regad)
    for (ii in 1:p2) {
      colnames(PCinterval1)[ii] <- paste(regad[ii, 1])
    }
    united = seq(1, 2 * p, 2)
    res <- matrix(0, p, 2)
    for (i in 1:p) {
      res[i, ] <- c(united[i], united[i] + 1)
    }
    zPC1 <- vector("list", p)
    for (i in 1:p) {
      zPC1[[i]] <- PCinterval1[, res[i, ]]
    }
    pc1 <- zPC1[[axes2[1]]]
    pc2 <- zPC1[[axes2[2]]]
    pc3 <- zPC1[[axes2[3]]]
    xmin1 <- min(pc1[, 1])
    xmin2 <- min(pc1[, 2])
    xxmin <- min(xmin1, xmin2)
    xmax1 <- max(pc1[, 1])
    xmax2 <- max(pc1[, 2])
    xxmax <- max(xmax1, xmax2)
    ymin1 <- min(pc2[, 1])
    ymin2 <- min(pc2[, 2])
    yymin <- min(ymin1, ymin2)
    ymax1 <- max(pc2[, 1])
    ymax2 <- max(pc2[, 2])
    yymax <- max(ymax1, ymax2)
    zmin1 <- min(pc3[, 1])
    zmin2 <- min(pc3[, 2])
    zzmin <- min(zmin1, zmin2)
    zmax1 <- max(pc3[, 1])
    zmax2 <- max(pc3[, 2])
    zzmax <- max(zmax1, zmax2)
    picto <- scatterplot3d(c(xxmin, xxmax), c(yymin, yymax), 
                           c(zzmin, zzmax), xlab = paste("PC", axes2[1], sep = "."), 
                           ylab = paste("PC", axes2[2], sep = "."), zlab = paste("PC", 
                                                                                 axes2[3], sep = "."), type = "n", main = "3D Factorial Plan", 
                           sub = nametable[k])
    cube <- rbind(c(1, 1, 1), c(2, 1, 1), c(2, 1, 2), c(1, 
                                                        1, 2), c(1, 1, 1), c(1, 2, 1), c(1, 2, 2), c(2, 2, 
                                                                                                     2), c(2, 2, 1), c(1, 2, 1), c(1, 2, 2), c(1, 1, 2), 
                  c(2, 1, 2), c(2, 2, 2), c(2, 2, 1), c(2, 1, 1))
    for (i in 1:n) {
      vec.x <- pc1[i, cube[, 1]]
      vec.y <- pc2[i, cube[, 2]]
      vec.z <- pc3[i, cube[, 3]]
      picto$points3d(vec.x, vec.y, vec.z, type = "l", lty = 1, 
                     col = i)
    }
    ltext <- Row.names
    textPoints <- cbind((pc1[, 1] + pc1[, 2])/2, (pc2[, 1] + 
                                                    pc2[, 2])/2, (pc3[, 1] + pc3[, 2])/2)
    text(picto$xyz.convert(textPoints), labels = ltext, col = 1:n)
    dev.new()
    xlimbis = range(tempGenPC)
    ylimbis = range(tempGenPC)
    n1 <- 1
    n2 <- n
    m1 <- n2 + 1
    m2 <- 2 * n2
    plot(tempGenPC[, axes[1]], tempGenPC[, axes[2]], xlab = lab.x, 
         ylab = lab.y, asp = 1, main = "Projection of table onto axes of the compromise.", 
         xlim = xlimbis, ylim = xlimbis, sub = nametable[k])
    text(tempGenPC[1:n2, axes[1]], tempGenPC[1:n2, axes[2]], 
         labels = Row.names, pos = 4, col = 1:n)
    abline(h = 0, col = "blue")
    abline(v = 0, col = "blue")
    rect(tempGenPC[n1:n2, axes[1]], tempGenPC[n1:n2, axes[2]], 
         tempGenPC[m1:m2, axes[1]], tempGenPC[m1:m2, axes[2]], 
         border = 1:n)
  }

  
  if(method=='hypercube')
  {
    k<-t
    n<-nrow(Variable[[1]])
    p<-length(Variable) 
    
    if(transformation==2)
    {
      for(j in 1:p) 
      {
        Variable[[j]] <- asin( sqrt(Variable[[j]]) )
      }
    }
    
    for(j in 1:p)
    {
      Variable[[j]]<-Centrage( (Variable[[j]]) )
    }
    if (is.null(score))# initialisation des score s'ils sont nuls
    {
      score<-vector('list',p)
      for (j in 1: p) 
        score[[j]]<-as.vector( c(  1:ncol( (Variable[[ j ]]) ) )  )
    }
    
    mean<-matrix(0,n,p)
    for(i in 1:n)
    {  for(j in 1:p)
      mean[i,j]<- sum( (( Variable[[j]])[i,] )*(score[[j]])  )
    }
    ##########################  Calcul des valeurs propres et vecteurs propres ################################
    
    V=(1/n)*t(mean)%*%mean
    colnames(V)<-paste("Variable ", 1:p, sep = '')
    rownames(V)<-paste("Variable ", 1:p, sep = '')
    #print('Matrice de variance covariance')
    #print(V)
    
    ###########################         DevStandard
    DevStandard<-matrix(0,n,p)
    for(i in 1:n)
    {  for(j in 1:p)
      DevStandard[i,j]<- sum( ( (( (Variable[[j]])[i,])^2 )*(score[[j]]) )  -(  sum( ( as.matrix(( Variable[[j]])[i,]  )*(score[[j]]) )  )^2 )/n )/(n-1)
    }
    Tmin<-mean-k*DevStandard
    Tmax<-mean+k*DevStandard
    
    auxil<-vector('list',p)
    for(j in 1:p)
      auxil[[j]]<-cbind(Tmin[,j],Tmax[,j])
    
    res<-auxil[[1]]
{
  for(s in 2:p)
    res<-cbind(res,auxil[[s]])             
}
    
    
    Vect<-svd(mean,nu=n, nv=p)$v
    Val<-svd(V,nu=p)$d
    val<-svd(V,nu=p)$d
    g=Variable[[1]]
    g<-as.data.frame(g)
    g<-Centrage(g)
    lo<-nrow(g)
    T=rbind(Tmin,Tmax)
    prepT=fintab(as.data.frame(T))
    gg=Hypercube(prepT)  
    gg<-as.data.frame(gg)
    gg<-Centrage(gg)
    glo<-nrow(gg)
    gdim<-ncol(gg)
    gn<-glo/(2^gdim)
    Mind<-cube(prepT)
    Compo<-Mind # intialisation
    for(i in 1:gn) 
    {
      Compo[[i]]<-as.matrix( Mind[[i]] )%*%Vect
    }
    
      
    CompoMin<-matrix(0,n,p)
    for(i in 1:n) 
      for(j in 1:p )
        CompoMin[i,j]<-min ( (Compo[[i]])[ ,j] )
    
    
    CompoMax<-matrix(0,n,p)
    for(i in 1:n) 
      for(j in 1:p )
        CompoMax[i,j]<-max( (Compo[[i]])[ ,j] )
    
   
        
    Pval<-Val/(sum(Val))
    
    Pval<-100*Pval
    Pval<-as.matrix(Pval)
    

	
	 cp1 <- Pval[axes[1]]
    cp2 <- Pval[axes[2]]
    
	
	 lab.x <- paste(" Axe n", axes[1], " (", cp1, "Percentage)",sep = "")
    lab.y <- paste(" Axe n", axes[2], " (", cp2, "Percentage)",sep = "")
    
    
    if(transformation==1)
    {            
    
      title <- "Factorial Plan"
     }
    
    if(transformation==2)
    {            
      title <- "Factorial Plan avec Transformation angulaire"
    }
    

    
    if( is.null(Row.names) )
    {
      Row.names<-rownames(Variable[[1]])
    }
    
    Row.namesInd=Row.names;
 
    if( is.null(xlim) &&  is.null(ylim) ) 
    {
      xlim<-range(cbind(CompoMin[,axes[1]], CompoMax[,axes[1]])) 
      ylim<-range(cbind(CompoMin[,axes[2]], CompoMax[,axes[2]])) 
  
      
    }
    
    
    if( is.null(xlegend) &&  is.null(ylegend) ) 
    {
      xlegend<-(min(   range(Variable[[axes[1]]][,1])   ) + min(   range(Variable[[axes[1]]][,1]))   )*0.75
      ylegend<-(min(   range(Variable[[axes[2]]][,2])  ) + min( range(Variable[[axes[2]]][,1]))   )*(0.15)
      
    }
    
    m1<-nrow(Variable[[1]])
    plot(CompoMax[,axes[1]],CompoMax[,axes[2]] , xlab = lab.x, ylab = lab.y, xlim=xlim , ylim=ylim , asp = 1, main = title) 
    text(CompoMax[,axes[1]],CompoMax[,axes[2]] ,labels = Row.names)
    abline(h=0,col='blue')
    abline(v=0,col='blue')
    points(CompoMin[,axes[1]],CompoMin[,axes[2]],pch=1:m1 )
    points(CompoMin[,axes[1]],CompoMin[,axes[2]])
    rect(CompoMin[,axes[1]],CompoMin[,axes[2]],CompoMax[,axes[1]],CompoMax[,axes[2]], border=1:m1 )
    ClassCompo<-mean%*%Vect
    Correl<-cor(mean,ClassCompo)
    
    if (is.null(Col.names))
    {
      Col.names<-paste("Variable ", 1:p, sep = '')
    }

	lab.x <- paste("Composante n", axes[1]," (", cp1, "Percentage)", sep = "")
    lab.y <- paste("Composante n", axes[2], " (", cp2, "Percentage)",sep = "")
    colnames(Correl)<-paste("Composante ", 1:p, sep = '')
    rownames(Correl)<-Col.names
    if (is.null(Col.names))
    {
      Col.names<-paste("Variable", 1:p, sep = '')
    }

    dev.new()
    
    title2 <- "Correlation circle."
    plot(Correl[,axes[1]],Correl[,axes[2]], xlab = lab.x, ylab = lab.y, xlim=c(-1.5,1.5) , ylim=c(-1.5,1.5), asp = 1, main = title2)
    text(Correl[,axes[1]],Correl[,axes[2]] ,labels =rownames(Correl),pos=4)
    segments(0,0,Correl[,axes[1]],Correl[,axes[2]] ,col=1:p )
    
    abline(h=0)
    abline(v=0)
    Pval<-Val/(sum(Val))
    Pval<-100*Pval
    Pval<-as.matrix(Pval)
    rownames(Pval)<-paste('Composante',1:p,'')

    colnames(Pval)<-c('Poucentage de variabilite des composantes ')
    colnames(Vect)<-paste("VecteurPropre ", 1:p, sep = '')
    CarreCompo<-as.matrix(ClassCompo)*as.matrix(ClassCompo)      
	 rownames(Pval) <- paste("Component n", 1:p, "")
        
		colnames(Pval) <- c("Percentage of variability")
        pval2 <- as.data.frame(matrix(NA, length(Val), 3))
        rownames(pval2) <- paste("comp", 1:length(Val))
        colnames(pval2) <- c("eigenvalue", "percentage of variance", 
            "cumulative percentage of variance")
        pval2[, "eigenvalue"] <- Val
        pval2[, "percentage of variance"] <- (Val/sum(Val)) * 
            100
        pval2[, "cumulative percentage of variance"] <- cumsum(pval2[, 
            "percentage of variance"])
    


    GenPC=rbind(CompoMin,CompoMax) 
    PCinterval0 <- fintab(as.data.frame(GenPC))
    regad <- ucad(p)
    p2 = nrow(regad)
    for (ii in 1:p2) {
      colnames(PCinterval0)[ii] <- paste(regad[ii, 1])
    }
    united = seq(1, 2 * p, 2)
    res <- matrix(0, p, 2)
    for (i in 1:p) {
      res[i, ] <- c(united[i], united[i] + 1)
    }
    zPC0 <- vector("list", p)
    for (i in 1:p) {
      zPC0[[i]] <- PCinterval0[, res[i, ]]
    }
    
    if (p>2)
    {
      dev.new()
      pc1 <- zPC0[[axes2[1]]]
      pc2 <- zPC0[[axes2[2]]]
      pc3 <- zPC0[[axes2[3]]]
      xmin1 <- min(pc1[, 1])
      xmin2 <- min(pc1[, 2])
      xxmin <- min(xmin1, xmin2)
      xmax1 <- max(pc1[, 1])
      xmax2 <- max(pc1[, 2])
      xxmax <- max(xmax1, xmax2)
      ymin1 <- min(pc2[, 1])
      ymin2 <- min(pc2[, 2])
      yymin <- min(ymin1, ymin2)
      ymax1 <- max(pc2[, 1])
      ymax2 <- max(pc2[, 2])
      yymax <- max(ymax1, ymax2)
      zmin1 <- min(pc3[, 1])
      zmin2 <- min(pc3[, 2])
      zzmin <- min(zmin1, zmin2)
      zmax1 <- max(pc3[, 1])
      zmax2 <- max(pc3[, 2])
      zzmax <- max(zmax1, zmax2)
      picto <- scatterplot3d(c(xxmin, xxmax), c(yymin, yymax), 
                             c(zzmin, zzmax), xlab = paste("PC", axes2[1], sep = "."), 
                             ylab = paste("PC", axes2[2], sep = "."), zlab = paste("PC", 
                                                                                   axes2[3], sep = "."), type = "n", main = "3D Factorial Plan")
      cube <- rbind(c(1, 1, 1), c(2, 1, 1), c(2, 1, 2), c(1, 
                                                          1, 2), c(1, 1, 1), c(1, 2, 1), c(1, 2, 2), c(2, 2, 
                                                                                                       2), c(2, 2, 1), c(1, 2, 1), c(1, 2, 2), c(1, 1, 2), 
                    c(2, 1, 2), c(2, 2, 2), c(2, 2, 1), c(2, 1, 1))
      for (i in 1:n) {
        vec.x <- pc1[i, cube[, 1]]
        vec.y <- pc2[i, cube[, 2]]
        vec.z <- pc3[i, cube[, 3]]
        
        picto$points3d(vec.x, vec.y, vec.z, type = "l", lty = 1, 
                       col = i)
      }
      ltext <- Row.names
      textPoints <- cbind((pc1[, 1] + pc1[, 2])/2, (pc2[, 1] + 
                                                      pc2[, 2])/2, (pc3[, 1] + pc3[, 2])/2)
      text(picto$xyz.convert(textPoints), labels = ltext, col = 1:n)
      
    }
    
    lxmin=CompoMin[,axes[1]]
    lxmax=CompoMax[,axes[1]]
    lymin=CompoMin[,axes[2]]
    lymax=CompoMin[,axes[2]]
    
      colnames(GenPC)<-paste('Prin', 1:p,sep='.')
      newGenPC=as.data.frame(GenPC[,axes])
    rownames(PCinterval0)=Row.names
 return(list(Correlation=Correl,VecteurPropre=Vect, Tableaumean=mean,PourCentageComposante=pval2,PCinterval=PCinterval0))

 }
  
if(method=='longueur')
{
  k<-t
  n<-nrow(Variable[[1]])
  p<-length(Variable)
  if(transformation==2)
  {
    for(j in 1:p) 
    {
      Variable[[j]] <- asin( sqrt(Variable[[j]]) )
    }
  }
  
  for(j in 1:p)
  {
    Variable[[j]]<-Centrage( (Variable[[j]]) )
  }
  
  if (is.null(score))
  {
    score<-vector('list',p)
    for (j in 1: p)
      
      score[[j]]<-as.vector( c(  1:ncol( (Variable[[ j ]]) ) )  )
  }
  
  mean<-matrix(0,n,p)
  for(i in 1:n)
  {  for(j in 1:p)
    mean[i,j]<- sum( (( Variable[[j]])[i,] )*(score[[j]])  )
  }

  V=(1/n)*t(mean)%*%mean
  colnames(V)<-paste("Variable ", 1:p, sep = '')
  rownames(V)<-paste("Variable ", 1:p, sep = '')
  DevStandard<-matrix(0,n,p)
  for(i in 1:n)
  {  for(j in 1:p)
    DevStandard[i,j]<- sum( ( (( (Variable[[j]])[i,])^2 )*(score[[j]]) )  -(  sum( ( as.matrix(( Variable[[j]])[i,]  )*(score[[j]]) )  )^2 )/n )/(n-1)
  }
  
  Tmin<-mean-k*DevStandard
  Tmax<-mean+k*DevStandard
  Tdistance<-(Tmax-Tmin)
  Vect<-svd(V,nv=p)$v
  Val<-svd(V,nu=p)$d
  Composante<-mean%*%Vect
  Correl<-cor(mean,Composante)
  Tdistance<-Centrage(Tdistance)
  Vdistance<-(1/n)*t(Tdistance)%*%Tdistance
  PTdistance<-Tdistance%*%Vect
  Pval<-Val/(sum(Val))
  Pval<-100*Pval
  Pval<-as.matrix(Pval)
   
   rownames(Pval) <- paste("component n", 1:p, "")
        
		colnames(Pval) <- c("Percentage of variability")
        pval2 <- as.data.frame(matrix(NA, length(Val), 3))
        rownames(pval2) <- paste("comp", 1:length(Val))
        colnames(pval2) <- c("eigenvalue", "percentage of variance", 
            "cumulative percentage of variance")
        pval2[, "eigenvalue"] <- Val
        pval2[, "percentage of variance"] <- (Val/sum(Val)) * 
            100
        pval2[, "cumulative percentage of variance"] <- cumsum(pval2[, 
            "percentage of variance"])


   cp1 <- Pval[axes[1]]
  cp2 <- Pval[axes[2]]
  lab.x <- paste(" Axe n", axes[1], " (", cp1, "Percentage)",sep = "")
  lab.y <- paste(" Axe n", axes[2], " (", cp2, "Percentage)",sep = "")
  title <- "Factorial Plan"
  if(transformation==2)
  {            
    title <- "Factorial Plan"
  }
  
  if( is.null(Row.names) )
  {
    Row.names<-rownames(Variable[[1]])
  }
  
  
  if( is.null(xlim) &&  is.null(ylim) ) 
  {
    xlim<-range( Composante[,axes[1]],PTdistance[,axes[1]] ) 
    ylim<-range(Composante[,axes[2]],PTdistance[,axes[2]] ) 
  }
  plot(Composante[,axes[1]],Composante[,axes[2]] , xlab = lab.x, ylab = lab.y, xlim=xlim , ylim=ylim , asp = 1, main = title, sub='Use of angular transformation') 
  text(Composante[,axes[1]],Composante[,axes[2]] ,labels = Row.names)
  abline(h=0,col='blue')
  abline(v=0,col='blue')
  n<-nrow(Variable[[1]])
  points(PTdistance[,axes[1]],PTdistance[,axes[2]])
  rect(Composante[,axes[1]],Composante[,axes[2]],PTdistance[,axes[1]],PTdistance[,axes[2]], border=1:n )
  if( is.null(xlegend) &&  is.null(ylegend) ) 
  {
    xlegend<-(min(   range(Variable[[axes[1]]][,1])   ) + min(   range(Variable[[axes[1]]][,1]))   )*0.75
    ylegend<-(min(   range(Variable[[axes[2]]][,2])  ) + min( range(Variable[[axes[2]]][,1]))   )*(0.15)
  }
  
  if (is.null(Col.names))
  {
    Col.names<-paste("Variable", 1:p, sep = '')
  }
  
   lab.x <- paste("Component n", axes[1], " (", cp1, "Percentage)",sep = "")
  lab.y <- paste("Component n", axes[2], " (", cp2, "Percentage)",sep = "")
  
  colnames(Correl)<-paste("Component n", 1:p, sep = '')
  rownames(Correl)<-Col.names
  
  if (is.null(Col.names))
  {
    Col.names<-paste("Variable ", 1:p, sep = '')
  }
  
  dev.new()
  title2 <- "HistPCA: Correlation circle."
  plot(Correl[,axes[1]],Correl[,axes[2]], xlab = lab.x, ylab = lab.y, xlim=c(-1.5,1.5) , ylim=c(-1.5,1.5), asp = 1, main = title2)
  text(Correl[,axes[1]],Correl[,axes[2]] ,labels =rownames(Correl),pos=4)
  segments(0,0,Correl[,axes[1]],Correl[,axes[2]] ,col=1:p )
  abline(h=0,col='red')
  abline(v=0,col='red')
  Pval<-Val/(sum(Val))
  Pval<-100*Pval
  Pval<-as.matrix(Pval)
  rownames(Pval)<-paste('Component n',1:p,'')
  colnames(Pval)<-c('Percentage of Component variability')
  colnames(Vect)<-paste("VecteurPropre n", 1:p, sep = '')
  
 
  colnames(mean)<-Col.names
    
if(proc==0)
{
  dev.new()
  par(mfrow=c(1,2))
  plot(Correl[,axes[1]],Correl[,axes[2]], xlab = lab.x, ylab = lab.y, xlim=c(-1.5,1.5) , ylim=c(-1.5,1.5), asp = 1, main = 'Correlation circle')
  text(Correl[,axes[1]],Correl[,axes[2]] ,labels =rownames(Correl),pos=4)
  segments(0,0,Correl[,axes[1]],Correl[,axes[2]] ,col=1:p )
  abline(h=0)
  abline(v=0)
 
    lab.x <- paste(" Axe n", axes[1], " (", cp1, "Percent)",sep = "")
  lab.y <- paste(" Axe n", axes[2], " (", cp2, "Percent)",sep = "")
  
  
  plot(Composante[,axes[1]],Composante[,axes[2]] , xlab = lab.x, ylab = lab.y, xlim=xlim , ylim=ylim , asp = 1, main = title) 
  text(Composante[,axes[1]],Composante[,axes[2]] ,labels = Row.names)
  abline(h=0,col='blue')
  abline(v=0,col='blue')
  
  n<-nrow(Variable[[1]])
  
  
  points(PTdistance[,axes[1]],PTdistance[,axes[2]])
  rect(Composante[,axes[1]],Composante[,axes[2]],PTdistance[,axes[1]],PTdistance[,axes[2]], border=1:n )
   GenPC=rbind(Composante,PTdistance)
   

   PCinterval0 <- fintab(as.data.frame(GenPC))
    regad <- ucad(p)
    p2 = nrow(regad)
    for (ii in 1:p2) {
      colnames(PCinterval0)[ii] <- paste(regad[ii, 1])
    }
    united = seq(1, 2 * p, 2)
    res <- matrix(0, p, 2)
    for (i in 1:p) {
      res[i, ] <- c(united[i], united[i] + 1)
    }
    zPC0 <- vector("list", p)
    for (i in 1:p) {
      zPC0[[i]] <- PCinterval0[, res[i, ]]
    }
   
  return(list(Correlation=Correl,VecteurPropre=Vect, Tableaumean=mean,PourCentageComposante=pval2,PCinterval=PCinterval0))

   }  
 

     else
  {
    
    dev.new()
	 llab.x <- paste("Axe n", axes[1], sep = "")
    llab.y <- paste("Axe n", axes[2], sep = "")
 

 dframe=as.data.frame(cbind(Composante,PTdistance))

     xlim=range(dframe[,axes[1]]);
     ylim=range(dframe[,axes[2]]);

    pro1<-GPA(dframe,scale=FALSE,group=c(p,p) )
    xlim<-range(pro1$Xfin)
    ylim<-range(pro1$Xfin)
    plot(pro1$Xfin[,axes,1], xlab=llab.x,ylab=llab.y,main='Factorial Plan With Procuste Analysis',xlim=xlim,ylim=ylim, pch='C',type='n')
    points(rbind(pro1$Xfin[,axes,1],pro1$Xfin[,axes,2]))
    text(pro1$Xfin[,axes,1] ,labels = Row.names)
    rect(pro1$Xfin[,axes[1],1],pro1$Xfin[,axes[2],1],pro1$Xfin[,axes[1],2], pro1$Xfin[,axes[2],2],border=1:n )
    abline(h=0)
    abline(v=0)
  GenPC=rbind(pro1$Xfin[,,1],pro1$Xfin[,,2])
   pp=ncol(GenPC);
   PCinterval0 <- fintab(as.data.frame(GenPC))
    regad <- ucad(pp)
     p2 = nrow(regad)
    for (ii in 1:p2) {
      colnames(PCinterval0)[ii] <- paste(regad[ii, 1])
    }
    united = seq(1, 2 * pp, 2)
    res <- matrix(0, pp, 2)
    for (i in 1:pp) {
      res[i, ] <- c(united[i], united[i] + 1)
    }
    zPC0 <- vector("list", pp)
    for (i in 1:pp) {
      zPC0[[i]] <- PCinterval0[, res[i, ]]
    }
     return(list(Correlation=Correl,VecteurPropre=Vect, Tableaumean=mean,PourCentageComposante=pval2,PCinterval=PCinterval0))
 }
  
  
}
 
}
