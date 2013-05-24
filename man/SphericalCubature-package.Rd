\name{SphericalCubature-package}
\alias{SphericalCubature-package}
\alias{SphericalCubature}
\docType{package}
\title{
Numerical integration over spheres and balls in n-dimensions; multivariate polar/spherical coordinates}
\description{
This package defines functions to integrate a (single) function f(x)=f(x[1],...,x[n])
over the unit sphere and balls in n-dimensional Euclidean space:
\deqn{\int_S f(s) ds \quad \quad \quad \quad \mathrm{and} \quad \quad \quad \int_B f(x) dx,}{(integral over S) f(s) ds,       (integral over B) f(x) dx,}
where the first integral is over the unit sphere S, an (n-1) dimensional surface,
and the second integral is over the unit ball B, an n dimensional solid.

There are three classes of methods:
\enumerate{
 \item exact methods for polynomials in any dimension (fast)
 \item a method due to Stroud for smooth integrands  on the sphere (in 
      dimensions n=3,4,...,16) (slower)
 \item adaptive methods for integrands with different behavior in different reasons (slowest) 
}
Methods 2 and 3 are approximations: like any numerical quadrature
algorithm, they may give bad results if the integrand changes abruptly on a small region.
This happens even in one dimension, and is more difficult to find and deal with
in higher dimensions.  (One attempt to handle is the 'split' versions of the adaptive methods,
functions adaptIntegrateSphereSplit and adaptIntegrateBallSplit,
where one can split the region of integration based on knowledge of the integrand.) 

An explicit goal was to get beyond the cases where n=2 or n=3, so some 
efficiency has been sacrificed. In all the methods, the higher the 
dimension n, the longer the compute time.  For methods 2 and 3, compute
times get noticeable when n > 5.  It was also necessary for the problem
that motivated this package to deal with integrands that had sharp
spikes, and that requires some sort of adaptive technique.

The package includes functions to convert to/from polar coordinates in higher dimensions.

This is a first attempt to provide methods for integrating
over spheres and balls in multiple dimensions.  One possible improvement is speed:
coding routines in C would give a significant increase in speed.  
Another possible extension is to include other multivariate integration
methods,  e.g. the package R2cuba.  This may provide a way to approximate 
higher dimensional integrals in some cases, if the integrand is well behaved.
Vector valued integrands are not supported; this could be incorporated in
the future.
%% Google search on "spherical cubature" finds more recent work

Constructive comments for improvements are welcome;  
actually implementing any suggestions will be dependent on time constraints.

Version history: 
\itemize{
\item 1.0.0 (2013-05-16) original package
\item 1.0.1 (2013-05-24) fix mistake in adaptIntegrateBallSplit, fix example in integratePolynomialSphere, add more documentation
}
}
\details{
\tabular{ll}{
Package: \tab SphericalCubature\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2013-05-24\cr
License: \tab GPL\cr
}
}
\author{John P. Nolan

Maintainer: John P. Nolan <jpnolan@american.edu>
}

\seealso{
\code{\link{integrateSpherePolynomial}}, \code{\link{integrateBallPolynomial}}, 
\code{\link{integrateSphereStroud11}}, \code{\link{sphereArea}}, \code{\link{ballVolume}},
\code{\link{polar2rect}}, \code{\link{rect2polar}}, 
\code{\link{adaptIntegrateSphere}}, \code{\link{adaptIntegrateSphereSplit}}, 
\code{\link{adaptIntegrateBall}}, \code{\link{adaptIntegrateBallSplit}}
}

\examples{
#  integral should just be the area of sphere in n dimensions
f1 <- function( x ) { return(1.0) }
n <- 3
sphereArea( n )
integrateSphereStroud11( f1, n )
p <- list(coef=1.0,k=matrix( rep(0L,n), nrow=1,ncol=n))
integrateSpherePolynomial( p )
adaptIntegrateSphere( f1, n )$value 

# test of polynomial integration
f2 <- function( x ) { return(x[1]^2) }
sphereArea(n)/n # exact answer
integrateSphereStroud11( f2, n )
p <- list(coef=1.0,k=matrix( c(2L,rep(0L,n-1)), nrow=1) )
integrateSpherePolynomial( p )
adaptIntegrateSphere( f2, n )$value
}
