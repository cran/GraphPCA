\name{Visu}
\alias{Visu}
%-
\title{
Visu
}
\description{
This function plots a scatterplot of histogram 
variables, using the R package ggplot2. 

}
\usage{
Visu(PC, Row.names = NULL, labs = NULL, axes = c(1, 2),size = 4)
}

\arguments{
  \item{PC}{
PCA data frame scores that contains   (PCxmin, PCxmax, PCymin, PCymax).
}
  \item{Row.names}{
Row names of concepts. }
  \item{labs}{
set the names of the axes.
}
  \item{axes}{
a length 2 vector specifying the components to plot.
}
 \item{size}{
Default value
}
}



\references{
Billard, L. and E. Diday (2006). Symbolic Data Analysis: conceptual statistics and data Mining.
Berlin: Wiley series in computational statistics.


Le-Rademacher J., Billard L. (2013). Principal component histograms from interval-valued
observations, Computational Statistics, v.28 n.5, p.2117-2138.

Makosso-Kallyth S. and Diday E. (2012). Adaptation of interval PCA to symbolic histogram 
variables, Advances in Data Analysis and Classification July, Volume 6, Issue 2, pp 147-159. 

}
\author{
Brahim Brahim Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and 
Sun Makosso-Kallyth <makosso.sun@gmail.com>
}


\examples{

data(movies)
ab = movies
ab = na.omit(ab)
Action = subset(ab,Action==1)
Action$genre = as.factor("Action")
Drama = subset(ab,Drama==1)
Drama$genre = as.factor("Drama")

Animation = subset(ab,Animation==1)
Animation$genre = as.factor("Animation")

Comedy = subset(ab,Comedy==1)
Comedy$genre = as.factor("Comedy")

Documentary = subset(ab,Documentary ==1)
Documentary $genre = as.factor("Documentary")


Romance = subset(ab,Romance==1)
Romance$genre = as.factor("Romance")

Short = subset(ab,Short==1)
Short$genre = as.factor("Short")

 ab = rbind(Action,Drama,Animation,Comedy,Documentary,Romance,Short)
 Hist1=PrepHistogram(X=sapply(ab[,3],unlist),Z=ab[,25],k=5)$Vhistogram
Hist2=PrepHistogram(X=sapply(ab[,4],unlist),Z=ab[,25],k=5)$Vhistogram
 Hist3=PrepHistogram(X=sapply(ab[,5],unlist),Z=ab[,25],k=5)$Vhistogram
Hist4=PrepHistogram(X=sapply(ab[,6],unlist),Z=ab[,25],k=5)$Vhistogram
 Hist5=PrepHistogram(X=sapply(ab[,7],unlist),Z=ab[,25],k=5)$Vhistogram
 
 ss1=Ridi(Hist1)$Ridit
 ss2=Ridi(Hist2)$Ridit
 ss3=Ridi(Hist3)$Ridit
 ss4=Ridi(Hist4)$Ridit
 ss5=Ridi(Hist5)$Ridit

 
HistPCA(list(Hist1,Hist2,Hist3,Hist4,Hist5),score=list(ss1,ss2,ss3,ss4,ss5))

res_pca=HistPCA(list(Hist1,Hist2,Hist3,Hist4,Hist5),score=list(ss1,ss2,ss3,ss4,ss5))
 
 Visu(res_pca$PCinterval)
}





