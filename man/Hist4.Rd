\name{Hist4}
\alias{Hist4}
\docType{data}
\title{
Horsepower variable.
}
\description{
 This  histogram datasets Hist4 come  from the aggregation of the fourth 
variable (horsepower)of the dataset collected 
by Ernesto Ramos and David Donoho (1982) which dealt with automobiles. 
The original data set  contains 406 observations and the following 8 variables. 
We create 13 new groups of observations (Year70,...,Year82) or clusters. After that,
for every cluster and all 6 other variables, we build histograms with 5 bins.
}
\usage{data(Hist4)}
\format{
    A data frame with the following histogram bins horsepower.1, horsepower.2, horsepower.3,
horsepower.4, horsepower.5.
  \describe{
    \item{\code{horsepower.1}}{a numeric vector}
    \item{\code{horsepower.2}}{a numeric vector}
    \item{\code{horsepower.3}}{a numeric vector}
    \item{\code{horsepower.4}}{a numeric vector}
    \item{\code{horsepower.5}}{a numeric vector}
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
data(Hist4)
}
