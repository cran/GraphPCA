\name{HistPCA}
\alias{HistPCA}
\title{HistPCA}
\description{
Performs a PCA of multiple tables of histogram variables.
}
\usage{
HistPCA(Variable = list, score = NULL, t = 1.1,
 axes = c(1, 2),Row.names = NULL, xlim = NULL,
 ylim = NULL, xlegend = NULL, ylegend = NULL,
 Col.names = NULL, transformation = 1,
 method = "hypercube", proc = 0,
 plot3d.table = NULL, axes2 = c(1, 2, 3))
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

}
\value{

 \item{Correlation}{Correlations between means of histogram and their principal components}
\item{Tablemean}{Table containing the average of histogram mean}
\item{eigenvector}{eigen vector of PCA of histogram mean}
\item{eigen_value_summary}{a matrix containing all the eigenvalues, 
the percentage of variance and the cumulative percentage of variance}
\item{PCinterval}{Data frame  containing the coordinates of 
the individuals on the principal axes}

}
\references{

Billard, L. and E. Diday (2006). Symbolic Data Analysis: conceptual statistics and data Mining.
Berlin: Wiley series in computational statistics.

Diday, E., Rodriguez O. and Winberg S. (2000). Generalization of the Principal Components Analysis
to Histogram Data, 4th European Conference on Principles and Practice of Knowledge Discovery in Data Bases,
 September 12-16, 2000, Lyon, France.

Donoho, D., & Ramos, E. (1982). Primdata: Data Sets for Use With PRIM-H. Version for
second (15-18, Aug, 1983) Exposition of Statistical
Graphics Technology, by American Statistical
Association.

Le-Rademacher J., Billard L. (2013). Principal component histograms from interval-valued
observations, Computational Statistics, v.28 n.5, p.2117-2138.

Makosso-Kallyth S. and Diday E. (2012). Adaptation of interval PCA to symbolic histogram 
variables, Advances in Data Analysis and Classification July, Volume 6, Issue 2, pp 147-159. 



}
\author{

Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and Sun Makosso-Kallyth <sun.makosso-kallyth@crchuq.ulaval.ca>

}


\examples{
data(Hist1)
data(Hist2)
data(Hist3)
data(Hist4)
data(Hist5)
data(Hist6)

example0=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'
),score=list(c(1:5),c(1:4),c(1:5),c(1:5),c(1:5),c(1:5)))
example0

example1=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'
))
example1


PC_example1=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'
))$PCinterval
example1

example2=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'),
transformation=2
)
example2


example3=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'),
method='longueur')
example3


example4=HistPCA(Variable=list
(Hist1,Hist2,Hist3,Hist4,Hist5
,Hist6),axes=c(1,2),Row.names=
paste('Year',70:82,sep='-'),
Col.names=c('mpg', 
'cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'),
method='longueur',proc=1)
example4

}