\name{SphericalCubature-package}
\alias{SphericalCubature-package}
\alias{SphericalCubature}
\docType{package}
\title{
Numerical integration over spheres and balls in n-dimensions; multivariate polar/spherical coordinates}
\description{Provides functions to integrate a function f(x)=f(x[1],...,x[n])
over the unit sphere and balls in n-dimensional Euclidean space:
\deqn{\int_S f(s) ds \quad \quad \quad \quad \mathrm{and} \quad \quad \quad \int_B f(x) dx,}{(integral over S) f(s) ds,       (integral over B) f(x) dx,}
where the first integral is over the unit sphere S, an (n-1) dimensional surface,
and the second integral is over the unit ball B, an n dimensional solid.

There are three classes of methods:
\enumerate{
 \item exact methods for polynomials in any dimension (fast)
 \item a method due to Stroud for smooth integrands  on the sphere (in 
      dimensions n=3,4,...,16) (slower)
 \item adaptive methods for integrands with different behavior in different regions (slowest) 
}
Methods 2 and 3 are approximations: like any numerical quadrature
algorithm, they may give inaccurate results if the integrand changes abruptly on a small region.
This happens even in one dimension, and is more difficult to find and deal with
in higher dimensions.  (One attempt to handle this difficulty is the 'split' versions of the adaptive methods,
functions \code{adaptIntegrateSpherePolarSplit} and \code{adaptIntegrateBallPolarSplit},
where one can split the region of integration based on knowledge of the integrand.) 

Version 1.1 of this package introduces new methods to integrate over spheres.  
Earlier versions used only polar coordinate representations of the sphere.  Now one
can use both polar representaions and triangulations of the sphere.  The latter has advantages in
some cases: it avoids the problems with polar coordinates giving regions that are sometimes
rectangles and sometimes triangles (which occurs at the poles), triangles can be approximately
equal area in any dimension, etc.  While adding these new routines,
names became confusing.  Apologies to anyone who has trouble because of this, but it seems
better in the long run to explicitly name functions based on their approach.  Hence
\code{adaptIntegrateSphere} has been renamed \code{adaptIntegrateSpherePolar} to indicate
that it uses polar coordinates, while the new function \code{adaptIntegrateSphereTri}
uses spherical triangles. 


An explicit goal was to get beyond the cases where n=2 or n=3, so some 
efficiency has been sacrificed. In all the methods, the higher the 
dimension n, the longer the compute time.  For methods 2 and 3, compute
times get noticeable when n > 5.  One application 
that motivated this package required the ability to work reliably with integrands that have sharp
spikes. That requires some sort of adaptive technique, with the possibility of telling
the integration algorithm where the spikes are.

The package also includes functions to convert to/from polar coordinates in higher dimensions.

This is a second attempt to provide methods for integrating
over spheres and balls in multiple dimensions.  One possible improvement is speed:
coding routines in C would give a significant increase in speed.  
Another possible extension is to include other multivariate integration
methods,  e.g. the package R2cuba.  This may provide a way to approximate 
higher dimensional integrals in some cases, if the integrand is well behaved.
Vector valued integrands are not supported in most functions (but are in \code{adaptIntegrateSphereTri}); 
this could be incorporated in the future.
%% Google search on "spherical cubature" finds more recent work

Please let me know if you find any mistakes. I will try to fix bugs promptly.

Constructive comments for improvements are welcome;  
actually implementing any suggestions will be dependent on time constraints.

Version history: 
\itemize{
\item 1.0.0 (2013-05-16) original package
\item 1.0.1 (2013-05-24) fix mistake in \code{adaptIntegrateBallPolarSplit}, fix example in \code{integratePolynomialSphere}, add more documentation
\item 1.0.2 (2013-12-18) add function \code{adaptIntegrateSphereTri3d} to integrate over spherical triangles in 3-dimensions
\item 1.1 (2016-05-14) add function \code{adaptIntegrateSphereTri} to integrate over hyperspherical triangles in n-dimensions. 
\item 1.2 (2016-07-23) improve \code{adaptIntegrateSphereTri} where the integration over octants worked, but integrals over other subdivisions did not.
     New version works over any subdivision that doesn't cross into different octants (this is checked). 
     Minor changes to documentation were made and more checks on input values were added.
\item 1.3 (2016-10-20)  Improve changes in version 1.2: remove the restricition on simplices in \code{adaptIntegrateSphereTri}.  
    The input simplices are analyzed and if a simplex is in more than one octant, it is automatically subdivided, giving a 
    list of simplices that exactly cover the same part of the sphere and respect each octant.
}
}

\author{John P. Nolan

Maintainer: John P. Nolan <jpnolan@american.edu>

This research was supported by an agreement with Cornell University, Operations 
Research & Information Engineering, under contract W911NF-12-1-0385 from the Army 
Research Development and Engineering Command.
}

\seealso{
\code{\link{integrateSpherePolynomial}}, \code{\link{integrateBallPolynomial}}, 
\code{\link{integrateSphereStroud11}}, \code{\link{sphereArea}}, \code{\link{ballVolume}},
\code{\link{polar2rect}}, \code{\link{rect2polar}}, 
\code{\link{adaptIntegrateSpherePolar}}, \code{\link{adaptIntegrateSpherePolarSplit}}, 
\code{\link{adaptIntegrateSphereTri}}, \code{\link{adaptIntegrateSphereTri3d}},
\code{\link{adaptIntegrateBallPolar}}, \code{\link{adaptIntegrateBallPolarSplit}}
}

\examples{
#  integral should just be the area of sphere in n dimensions
f1 <- function( x ) { return(1.0) }
n <- 3
sphereArea( n )
integrateSphereStroud11( f1, n )
p <- list(coef=1.0,k=matrix( rep(0L,n), nrow=1,ncol=n))
integrateSpherePolynomial( p )
adaptIntegrateSpherePolar( f1, n )$integral 

# test of polynomial integration
f2 <- function( x ) { return(x[1]^2) }
sphereArea(n)/n # exact answer
integrateSphereStroud11( f2, n )
p <- list(coef=1.0,k=matrix( c(2L,rep(0L,n-1)), nrow=1) )
integrateSpherePolynomial( p )
adaptIntegrateSpherePolar( f2, n )$integral

}
