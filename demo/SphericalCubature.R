# more examples for SphericalCubature, some are slow
######################################################################
# different dimensions
f <- function( s ) { return( s[1]^2 ) }
for (n in 2:8) {
  b <- adaptIntegrateSphereTri( f, n, Orthants(n,positive.only=TRUE) )
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
  a <- adaptIntegrateSphereTri( f1, n, maxEvals=1000000 )
  cat("  ",a$integral)
  if(n==3) {
    b <- adaptIntegrateSphereTri3d( f1, Orthants(n) )
    cat("  ",b$integral)
  }
  cat('\n')
}

# check that different subdivisions of sphere work correctly
for (n in 2:3) { 
  for (k in 0:3) {
    s1 <-  UnitSphere( n, k )
    sphere <- aperm( s1$S, c(2,1,3) )
    a <- adaptIntegrateSphereTri( f1, n, partitionInfo=TRUE )
    cat(n,k,a$integral,sphereArea( n ),"\n") # (calculated value, exact answer)
  }
}

# test of polynomial integration f(s)=s[1]^2
f5 <- function( s ) { return( s[1]^2 ) }
n <- 3
S <- Orthants( n, positive.only=TRUE )
a <- adaptIntegrateSphereTri( f5, n, S )
a$integral
sphereArea(n)/(2^n*n) - a$integral  # error=exact-cubature approximation

# increase accuracy by specifying a smaller tolerance
a <- adaptIntegrateSphereTri( f5, n, S, tol=1.0e-12 ) 
a$message # but see that we exceed maximum number of evaluations
sphereArea(n)/(2^n*n) - a$integral  # error improved, but not what we asked for
a <- adaptIntegrateSphereTri( f5, n, S, maxEvals=200000, tol=1.0e-12 ) # increase accuracy
a$message
sphereArea(n)/(2^n*n) - a$integral  # still too many evals, but better accuracy

# integrate over a subdivision of the positve orthant; note that this improves accuracy 
#       and decreases the number of function evaluations
S <- aperm( UnitSphere(n,k=2,positive.only=TRUE)$S, c(2,1,3) )
a <- adaptIntegrateSphereTri( f5, n, S) # increase accuracty
sphereArea(n)/(2^n*n) - a$integral  # error=exact-cubature approximation


# check that simplices crossing orthants gives correct answer
b <- sqrt(2)/2
badS <- matrix( c(b,b,  -b, b), nrow=2, ncol=2) # 2d example
adaptIntegrateSphereTri( f1, n=2, badS ) 

