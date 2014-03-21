\name{Hist2}
\alias{Hist2}
\docType{data}
\title{
Number cylinders
}
\description{
  This  histogram datasets Hist2 comes  from the aggregation of the second 
variable (number cylinders)of the dataset collected 
by Ernesto Ramos and David Donoho (1982) which dealt with automobiles. 
The original data set  contains 406 observations and the following 8 variables. 
We create 13 new groups of observations (Year70,...,Year82) or clusters. After that,
for every cluster and all 6 other variables, we build histograms with 4 bins.

}
\usage{data(Hist2)}
\format{
Data frame with the following histogram bins cylinders.1, cylinders.2, cylinders.3, cylinders.4.
  \describe{
    \item{\code{cylinders.1}}{a numeric vector}
    \item{\code{cylinders.2}}{a numeric vector}
    \item{\code{cylinders.3}}{a numeric vector}
    \item{\code{cylinders.4}}{a numeric vector}
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
data(Hist2)
}

