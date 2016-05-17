# more examples for SphericalCubature, some are slow
######################################################################
# different dimensions
f <- function( s ) { return( s[1]^2 ) }
for (n in 2:8) {
  b <- adaptIntegrateSphereTri( f, Octants(n,positive.only=TRUE) )
  cat( n, b$integral, sphereArea(n)/(n*2^n),"  ratio=",
    b$integral/( sphereArea(n)/(n*2^n)),"\n")
}
######################################################################
# spherical cubature - evaluate multiple ways
f1 <- function( x ) { return(1.0) }
for (n in 2:4) {
  if(n==2) cat("n   exact     stroud     poly      polar        Tri       Tri3d\n")
  cat(n,sphereArea( n ))
  if (n > 2) {cat("  ",integrateSphereStroud11( f1, n ))} else {cat("           ")}
  p <- list(coef=1.0,k=matrix( rep(0L,n), nrow=1,ncol=n))
  cat("  ",integrateSpherePolynomial( p ))
  cat("  ",adaptIntegrateSpherePolar( f1, n )$integral)
  a <- adaptIntegrateSphereTri( f1, Octants(n), maxEvals=1000000 )
  cat("  ",a$integral)
  if(n==3) {
    b <- adaptIntegrateSphereTri3d( f1, Octants(n) )
    cat("  ",b$integral)
  }
  cat('\n')
}
