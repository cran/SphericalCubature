\name{SphericalCubature-package}
\alias{SphericalCubature-package}
\alias{SphericalCubature}
\docType{package}
\title{
Numerical integration over spheres and balls in n-dimensions; multivariate polar/spherical coordinates}
\description{Provides functions to integrate a function f(x)=f(x[1],...,x[n])
over the unit sphere and unit ball in n-dimensional Euclidean space:
\deqn{\int_S f(s) ds \quad \quad \quad \quad \mathrm{and} \quad \quad \quad \int_B f(x) dx,}{(integral over S) f(s) ds,       (integral over B) f(x) dx,}
where the first integral is over the unit sphere S, an (n-1) dimensional surface,
and the second integral is over the unit ball B, an n dimensional solid.
Vector valued integrands are allowed in some functions; see functions with an argument named \code{fDim}.

Integration over general sphers and balls can be done by adjusting the integrand function: if S2 is sphere with center x0 and radius R and B2 is a ball with
same center and radius, a simple change of variable shows that 
\deqn{\int_{S2} f(t) dt = \int_S R^{n-1} f(x0+R s) ds  \quad \quad \quad \mathrm{and}  \quad \quad \int_{B2} f(y) dy = \int_B R^n f(x0+R x) dx.}{
(integral over S2) f(t) dt = (integral over S) R^{n-1} * f(x0+R*s) ds,      (integral over B2) f(y) dy = (integral over B) R^n * f(x0+R*x) dx.}
See the example in \code{\link{adaptIntegrateBallTri}}.

The package also includes functions to convert to/from polar coordinates in higher dimensions.

There are three cubature methods:
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
functions \code{adaptIntegrateSpherePolarSplit()} and \code{adaptIntegrateBallPolarSplit()},
where one can split the region of integration based on knowledge of the integrand.
This can also be done in functions \code{adaptIntegrateSphereTri()} and \code{adaptIntegrateBallTri()}
by specifying an appropriate partition in the simplices.)

It is expected that these methods will yield several significant digits, but not
many digits of precision.  This seems to be the state of the art in multivariate integration.

Version 1.1 of this package introduces new methods to integrate over spheres.  
Earlier versions used only polar coordinate representations of the sphere.  Now one
can use both polar representaions and triangulations of the sphere.  The latter has advantages in
some cases: it avoids the problems with polar coordinates giving regions that are sometimes
rectangles and sometimes triangles (which occurs at the poles), triangles can be approximately
equal area in any dimension, etc.  While adding these new routines,
names became confusing.  Apologies to anyone who has trouble because of this, but it seems
better in the long run to explicitly name functions based on their approach.  Hence
\code{adaptIntegrateSphere()} has been renamed \code{adaptIntegrateSpherePolar()} to indicate
that it uses polar coordinates, while the functions \code{adaptIntegrateSphereTri()} and 
\code{adaptIntegrateBallTri()} uses spherical triangles. 


An explicit goal was to get beyond the cases where n=2, so some 
efficiency has been sacrificed. In all the methods, the higher the 
dimension n, the longer the compute time.  For methods 2 and 3, compute
times get noticeable when n > 5.  One application 
that motivated this package required the ability to work reliably with integrands that have 
spikes. That requires some sort of adaptive technique, with the possibility of telling
the integration algorithm where the spikes are.


This package is an attempt to provide methods for integrating
over spheres and balls in multiple dimensions, not a final answer.  One possible improvement is speed:
coding routines in C would give a significant increase in speed.  
Another possible extension is to include other multivariate integration
methods,  e.g. the package R2cuba.  This may provide a way to approximate 
higher dimensional integrals in some cases, if the integrand is well behaved.
%% Google search on "spherical cubature" finds more recent work

Please let me know if you find any mistakes. Constructive comments for improvements are welcome.  
Fixing bugs or implementing suggestions will be dependent on my workload.

Version history: 
\itemize{
\item 1.0.0 (2013-05-16) original package
\item 1.0.1 (2013-05-24) fix mistake in \code{adaptIntegrateBallPolarSplit}, fix example in \code{integratePolynomialSphere}, add more documentation
\item 1.0.2 (2013-12-18) add function \code{adaptIntegrateSphereTri3d} to integrate over spherical triangles in 3-dimensions
\item 1.1 (2016-05-14) add function \code{adaptIntegrateSphereTri} to integrate over spherical triangles in n-dimensions. 
\item 1.2 (2016-07-23) improve \code{adaptIntegrateSphereTri} where the integration over octants worked, but integrals over other subdivisions did not.
     New version works over any subdivision that doesn't cross into different octants (this is checked). 
     Minor changes to documentation were made and more checks on input values were added.
\item 1.3 (2017-09-16)  Improve changes in version 1.2: remove the restricition on simplices in \code{adaptIntegrateSphereTri}:  
    the input simplices are analyzed and if a simplex is in more than one orthant, it is automatically subdivided, giving a 
    list of simplices that exactly cover the same part of the sphere and respect each orthant.
    Fix adaptIntegrateSphericalCubatureTri to correctly pass optional arguments to the integrand function. 
    Change the word "octant" to "orthant" throughout the code to stress that the code works in higher dimensions.
\item 1.4 (2017-09-16) minor tweaks
\item 1.5 (2021-01-04) (a) Fix bug in function \code{rect2polar()} (conversion from rectangular to 
    polar coordinates) that occurred along an axis.  Thanks to Prof. Eckard Liebscher for 
    finding this mistake. (b) Add function \code{adaptIntegrateBallSimplices()}.  
    (c) Use the new function \code{adaptIntegrateVectorFunc()} from package
    SimplicialCubature to handle the 2-dimensional case of integrating over a sphere. 
    (d) Add a required argument to function \code{adaptIntegrateSphereTri()}.
    (e) Add function \code{adaptIntegrateBallRadial()} to integrate a radial function over the unit ball. 
    (f) Recode function \code{Orthants()} to be a standalone
    function, not requiring package \code{mvmesh}. (g) add to documentation.
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
\code{\link{adaptIntegrateBallPolar}}, \code{\link{adaptIntegrateBallPolarSplit}},
\code{\link{adaptIntegrateBallTri}}, 
\code{\link{adaptIntegrateBallRadial}}
}

\examples{
#  integral should just be the area of sphere in n dimensions
f1 <- function( x ) { return(1.0) }
n <- 3
adaptIntegrateBallTri( f1, n )   # exact answer = volume of ball = (4/3)*pi = 4.18879

# other methods
integrateSphereStroud11( f1, n )
adaptIntegrateSpherePolar( f1, n )$integral 
# exact value for a polynomial
p <- list(coef=1.0,k=matrix( rep(0L,n), nrow=1,ncol=n))
integrateSpherePolynomial( p )

# test of exact polynomial integration
f2 <- function( x ) { return(x[1]^2) }
sphereArea(n)/n # exact answer
integrateSphereStroud11( f2, n )
p <- list(coef=1.0,k=matrix( c(2L,rep(0L,n-1)), nrow=1) )
integrateSpherePolynomial( p )
adaptIntegrateSpherePolar( f2, n )$integral

# integration over a ball
adaptIntegrateSphereTri( f1, n ) # exact answer = surface area of sphere = 4*pi = 12.56637

# vector valued integrands
f3 <- function( x ) { c( x[1]^2, x[2]^3 ) }
adaptIntegrateBallTri( f3, n, fDim=2 )
adaptIntegrateSphereTri( f3, n, fDim=2 )

# radial function
g <- function( x ) { sum(x^2) }
adaptIntegrateBallRadial( g, n=3 )

# for more examples enter:   demo(SphericalCubature)
}
