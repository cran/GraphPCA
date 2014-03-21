\name{Hist1}
\alias{Hist1}
\docType{data}
\title{
Miles per gallon
}
\description{
This  histogram datasets Hist1 comes  from the aggregation of the first
variable (Miles per gallon)of the dataset collected 
by Ernesto Ramos and David Donoho (1982) which dealt with automobiles. 
The original data set  contains 406 observations and the following 8 variables. 
We create 13 new groups of observations (Year70,...,Year82) or clusters. After that,
for every cluster and all 6 other variables, we build histograms with 5 bins Every histogram has 5 bins.
}
\usage{data(Hist1)}
\format{
Data frame with the following histogram bins mpg.1, mpg.2, mpg.3, mpg.4, mpg.5. 
  \describe{
    \item{\code{mpg.1}}{a numeric vector}
    \item{\code{mpg.2}}{a numeric vector}
    \item{\code{mpg.3}}{a numeric vector}
    \item{\code{mpg.4}}{a numeric vector}
    \item{\code{mpg.5}}{a numeric vector}
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
data(Hist1)
}

