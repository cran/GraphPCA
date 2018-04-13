\name{HistPCA}
\alias{HistPCA}
\title{HistPCA}
\description{
Performs a PCA of multiple tables of histogram variables.
}
\usage{
HistPCA(Variable = list, score = NULL, t = 1.1, axes = c(1, 2), 
        Row.names = NULL, xlim = NULL, ylim = NULL, xlegend = NULL, ylegend = NULL,
        Col.names = NULL, transformation = 1, method = "hypercube", proc = 0,
        plot3d.table = NULL, axes2 = c(1, 2, 3), ggplot = 1)
}

\arguments{
  \item{Variable}{
List of all data frames containing initial histogram variable. Every histogram is a data frames
and every columns of data frame contains histogram bins.  
}
  \item{score}{
List of bins score of every histogram variable.
By default these scores are the ranks of histogram bins.
}
  \item{t}{
t is a real number used for transforming histogram to interval via Tchebytchev's inequality. 
By default, t=1.1.
}
  \item{axes}{
a length 2 vector specifying the components to plot }

  \item{Row.names}{
Retrieve or set the row names of a matrix-like object.
}
  \item{xlim}{
range for the plotted "x" values, defaulting to the range of the finite values of "x".
}
  \item{ylim}{
range for the plotted "y" values, defaulting to the range of the finite values of "y".
}
  \item{xlegend}{
This function could be used to add legends to plots.
}
  \item{ylegend}{
This function could be used to add legends to plots.
}
  \item{Col.names}{
Retrieve or set the row names of a matrix-like object.
}
  \item{transformation}{
type of tranformation for data. If transformation=2, angular is used.
}
  \item{method}{
method used (method='hypercube',method='longueur')
}
  \item{proc}{
option valid when method='longueur'. If proc=1, the procuste analysis is used.
}
  \item{plot3d.table}{
specification for the scatterplot3d. if plot3d.table=1, the scatterplot3d will
appear.
}
  \item{axes2}{
a length 2 vector specifying the components to plot
}
  \item{ggplot}{
%%     ~~Describe \code{ggplot} here~~
}
}
\details{
See Examples
}
\value{

  \item{Correlation}{Correlations between means of histogram and their principal components}
 \item{Tableaumean}{Table containing the average of histogram mean}
 \item{VecteurPropre}{eigen vector of PCA of histogram mean}
\item{PourCentageComposante}{a matrix containing all the eigenvalues, 
the percentage of variance and the cumulative percentage of variance}
\item{PCinterval}{Data frame  containing the coordinates of 
the individuals on the principal axes}

}
\references{

Billard, L. and E. Diday (2006). Symbolic Data Analysis: conceptual statistics and data Mining.
Berlin: Wiley series in computational statistics.

Diday, E., Rodriguez O. and Winberg S. (2000). Generalization of the Principal Components Analysis
to Histogram Data, 4th European Conference on Principles and Practice of Knowledge Discovery in Data Bases, September 12-16, 2000, Lyon, France.

Donoho, D.,  Ramos, E. (1982). Primdata: Data Sets for Use With PRIM-H. Version for
second (15-18, Aug, 1983) Exposition of Statistical
Graphics Technology, by American Statistical
Association.

Le-Rademacher J., Billard L. (2013). Principal component histograms from interval-valued
observations, Computational Statistics, v.28 n.5, p.2117-2138.

Makosso-Kallyth S. and Diday E. (2012). Adaptation of interval PCA to symbolic histogram 
variables, Advances in Data Analysis and Classification July, Volume 6, Issue 2, pp 147-159. 



}
\author{

Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and Sun Makosso-Kallyth <makosso.sun@gmail.com>

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



