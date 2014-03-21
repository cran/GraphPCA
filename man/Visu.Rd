\name{Visu}
\alias{Visu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visu
}
\description{
This function plots a scatterplot of histogram 
variables, using the R package ggplot. 

}
\usage{
Visu(PC, Row.names = NULL, labs = NULL, axes = c(1, 2))
}

\arguments{
  \item{PC}{
data frame with for columns (xmin, xmax, ymin, ymax) containing symbolic 2 
symbolic interval variables.
}
  \item{Row.names}{
Retrieve or set the row names of a matrix-like object. }
  \item{labs}{
set the names of the axes.
}
  \item{axes}{
a length 2 vector specifying the components to plot.
}
}



\references{
Billard, L. and E. Diday (2006). Symbolic Data Analysis: conceptual statistics and data Mining.
Berlin: Wiley series in computational statistics.

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
Brahim Brahim Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and 
Sun Makosso-Kallyth <sun.makosso-kallyth@crchuq.ulaval.ca>
}


\examples{
data(Hist1)
data(Hist2)
data(Hist3)
data(Hist4)
data(Hist5)
data(Hist6)


PC_example1=HistPCA(Variable=
list(Hist1,Hist2,Hist3,Hist4,
Hist5,Hist6),axes=c(1,2),
Row.names=paste('Year',
70:82,sep='-'),Col.names=
c('mpg','cylinders', 
'displacement',
'horsepower',
'weight',
'acceleration'
))$PCinterval


Visu(PC_example1, axes=c(1,2),
Row.names=rownames(PC_example1))

}
