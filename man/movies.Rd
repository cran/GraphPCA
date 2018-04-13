\name{movies}
\alias{movies}
\docType{data}
\title{movies
}
\description{
A movies data frame with
}
\usage{data("movies")}
\format{
  A data frame with 58788 observations on the following 24 variables.
  \describe{
    \item{\code{title}}{a character vector}
    \item{\code{year}}{a numeric vector}
    \item{\code{length}}{a numeric vector}
    \item{\code{budget}}{a numeric vector}
    \item{\code{rating}}{a numeric vector}
    \item{\code{votes}}{a numeric vector}
    \item{\code{r1}}{a numeric vector}
    \item{\code{r2}}{a numeric vector}
    \item{\code{r3}}{a numeric vector}
    \item{\code{r4}}{a numeric vector}
    \item{\code{r5}}{a numeric vector}
    \item{\code{r6}}{a numeric vector}
    \item{\code{r7}}{a numeric vector}
    \item{\code{r8}}{a numeric vector}
    \item{\code{r9}}{a numeric vector}
    \item{\code{r10}}{a numeric vector}
    \item{\code{mpaa}}{a character vector}
    \item{\code{Action}}{a numeric vector}
    \item{\code{Animation}}{a numeric vector}
    \item{\code{Comedy}}{a numeric vector}
    \item{\code{Drama}}{a numeric vector}
    \item{\code{Documentary}}{a numeric vector}
    \item{\code{Romance}}{a numeric vector}
    \item{\code{Short}}{a numeric vector}
  }
}
\details{
Initial movies data frame on which Histogram variables are built/ 
}
\source{
https://cran.r-project.org/web/packages/ggplot2movies/index.html
}
\references{
Makosso-Kallyth, Sun; Diday, Edwin. Adaptation of interval PCA to symbolic histogram variables. 
Advances in Data Analysis and Classification. Volume 6. n 2. 2012. pages 147-159. Springer. 
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
\keyword{datasets}
