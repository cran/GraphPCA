\name{PrepHistogram}
\alias{PrepHistogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function transform standard variables into histogram-valued variable (PrepHistogram)
}
\description{
This function transforms quantitative variable into histogram-valued variable.
}
\usage{
PrepHistogram(X, Z = NULL, k = 3,group=NULL)
}
\arguments{
  \item{X}{
 \code{X}  quantitative variable that need to be transformed into symbolic histogram 

}
  \item{Z}{
 \code{Z} categorical variable that need to be used for  the purpose of clustering

}
  \item{k}{
 \code{k} number of bins of histogram-valued variables
}

 \item{group}{
  Data group
}

}

\value{
Returns a list including : 
  \item{mido }{return the class centres of class and the min and the max of histograms}
  \item{Vhistogram}{dataframe containing the relative frequency of histogram variables}

}
\references{
Makosso-Kallyth, Sun; Diday, Edwin. Adaptation of interval PCA to symbolic histogram variables. 
Advances in Data Analysis and Classification. Volume 6. n 2. 2012. pages 147-159. Springer. 
}
\author{
Brahim Brahim Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and 
Sun Makosso-Kallyth <makosso.sun@gmail.com>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
http://www.ivisualizations.com
}
\examples{

###  example1 from iris data

## preparation of histogram-valued variables (k= 3 bins)

Sepal.LengthHistogram=PrepHistogram(X=iris[,1],Z=iris[,5])$Vhistogram
Sepal.WidthHistogram=PrepHistogram(X=iris[,2],Z=iris[,5])$Vhistogram
Petal.LengthHistogram=PrepHistogram(X=iris[,3],Z=iris[,5])$Vhistogram
Petal.WidthHistogram=PrepHistogram(X=iris[,4],Z=iris[,5])$Vhistogram

############################     Hitsogram PCA   #################################
HistPCA(Variable=list(Sepal.LengthHistogram,Sepal.WidthHistogram,
        Petal.LengthHistogram,Petal.WidthHistogram),
        Row.names=names(table(iris[,5])),
		Col.names=colnames(iris)[1:4])


###  example2 from iris data

## preparation of histogram-valued variables (k= 4 bins)

Sepal.LengthHistogram=PrepHistogram(X=iris[,1],Z=iris[,5],k=2)$Vhistogram
Sepal.WidthHistogram=PrepHistogram(X=iris[,2],Z=iris[,5],k=2)$Vhistogram
Petal.LengthHistogram=PrepHistogram(X=iris[,3],Z=iris[,5],k=2)$Vhistogram
Petal.WidthHistogram=PrepHistogram(X=iris[,4],Z=iris[,5],k=2)$Vhistogram


############################     Hitsogram PCA   #################################

HistPCA(Variable=list(Sepal.LengthHistogram,Sepal.WidthHistogram,
Petal.LengthHistogram,Petal.WidthHistogram),
Row.names=names(table(iris[,5])),Col.names=colnames(iris)[1:4])
}
