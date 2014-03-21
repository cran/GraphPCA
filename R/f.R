f <-
function(dimension) {
   rorthogonal <- function(n = 2, dimension) {
       # generate `n` random orthonormal vectors
      # in `dimension` dimensional space
      x <- matrix(rnorm(n*dimension), nrow=n, ncol=dimension)
       t(eigen(crossprod(x))$vectors)
  }
  cube <- as.matrix(do.call(expand.grid, rep(list(c(-1,1)),dimension)))
  # the points pairs with distance 2 are connected by edges
  ij <- which(abs(as.matrix(dist(cube))-2)<1e-7, arr.ind=TRUE)
  # only need edge in one direction
  ij <- ij[ij[,1]<ij[,2],]
  basis <- rorthogonal(n=2, dim=dimension)
   # p is projection of vectices to subspace of basis
   p <- t(lm.fit(x=t(basis), y=t(cube))$coef)
   # plot vertices
   plot(p[,1],p[,2],asp=1)
  # plot edges
   segments(p[,1][ij[,1]], p[,2][ij[,1]], p[,1][ij[,2]], p[,2][ij[,2]])
return(list(p=p,ij=ij,cube=cube,basis=basis))
 }
