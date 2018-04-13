\name{Ridi3}
\alias{Ridi3}
%- Mantel Hansen Scores
\title{
Normalized Mantel Hansen's Scores computation using cumulative distribution function
}
\description{
 This function computes Normalized Mantel Hansen's scores
}
\usage{
Ridi3(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
 histogram-valued variable built via PrepHistogram function.
}
}
\details{
This function computes Normalized Mantel Hansen's scores
}
\value{
Returns a list including : 
  \item{p }{bins number of histograms}
  \item{entree }{Means table of Mantel's Hansen scores}
  \item{Ridit}{Mantel's Hansen scores}
}
\references{
Makosso-Kallyth, Sun; Diday, Edwin. Adaptation of interval PCA to symbolic histogram variables. 
Advances in Data Analysis and Classification. Volume 6. n 2. 2012. pages 147-159. Springer.
}
\author{
Brahim Brahim Brahim Brahim <brahim.brahim@bigdatavisualizations.com> and 
Sun Makosso-Kallyth <makosso.sun@gmail.com>
}
\note{
Computes Normalized Mantel Hansen's scores
}



\examples{

Hist1=PrepHistogram(X=iris[,1],Z=iris[,5])$Vhistogram
Hist2=PrepHistogram(X=iris[,2],Z=iris[,5])$Vhistogram
Hist3=PrepHistogram(X=iris[,3],Z=iris[,5])$Vhistogram
Hist4=PrepHistogram(X=iris[,4],Z=iris[,5])$Vhistogram

sss1=Ridi3(Hist1)$Ridit
sss2=Ridi3(Hist2)$Ridit
sss3=Ridi3(Hist3)$Ridit
sss4=Ridi3(Hist4)$Ridit

}

% Add one or more standard keywords, see file 'KEYWORDS' in the

