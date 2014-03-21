\name{Hist3}
\alias{Hist3}
\docType{data}
\title{
Engine displacement
}
\description{
  This  histogram datasets Hist3 comes  from the aggregation of the third 
variable (engine displacement (cu. inches))of the dataset collected 
by Ernesto Ramos and David Donoho (1982) which dealt with automobiles. 
The original data set  contains 406 observations and the following 8 variables. 
We create 13 new groups of observations (Year70,...,Year82) or clusters. After that,
for every cluster and all 6 other variables, we build histograms with 5 bins.
}
\usage{data(Hist3)}
\format{
  A data frame with the following histogram bins displacement.1, displacement.2,
 displacement.3, displacement.4 and displacement.5.
  \describe{
    \item{\code{displacement.1}}{a numeric vector}
    \item{\code{displacement.2}}{a numeric vector}
    \item{\code{displacement.3}}{a numeric vector}
    \item{\code{displacement.4}}{a numeric vector}
    \item{\code{displacement.5}}{a numeric vector}
  }
}

\source{
http://lib.stat.cmu.edu/datasets/cars.data
}
\references{
Donoho, D., & Ramos, E. (1982). Primdata: Data
Sets for Use With PRIM-H, (1982). Version for
second (15-18, Aug, 1983) Exposition of Statistical
Graphics Technology, by American Statistical
Association.

}
\examples{
data(Hist3)
}

