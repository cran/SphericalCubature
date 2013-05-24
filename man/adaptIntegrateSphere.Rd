\name{adaptIntegrateSphere}
\alias{adaptIntegrateSphere}
\alias{adaptIntegrateSphereSplit}
\alias{adaptIntegrateBall}
\alias{adaptIntegrateBallSplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adaptive integration over sphere or ball in n-dimensions
}
\description{
Approximate the integral over the sphere or ball in n-dimensions.  Can also
integrate over sectors of the sphere/ball, see details.  These functions
will be slow, but may be necessary to get accurate answers if the integrand
function f(x) is not well-behaved.  If the integrand changes rapidly
in certain regions, the basic routines adaptIntegrateSphere and 
adaptIntegrateBall will likely miss these abrupt changes and give
inaccurate results.  For cases where the location of the rapid changes are known,
the functions adaptIntegrateSphereSplit and
adaptIntegrateBallSplit allow you to split the region of integration
and capture those changes.
}
\usage{
adaptIntegrateSphere(f, n, lowerLimit = rep(0, n - 1), 
    upperLimit = c(rep(pi, n - 2), 2 * pi), tol = 1e-05, ...) 
adaptIntegrateSphereSplit(f, n, xstar, width = 0, lowerLimit = rep(0, n - 1), 
    upperLimit = c(rep(pi, n - 2), 2 * pi), tol = 1e-05, ...)

adaptIntegrateBall(f, n, lowerLimit = rep(0, n - 1), 
    upperLimit = c(rep(pi, n - 2), 2 * pi), R = c(0, 1), tol = 1e-05, ...)
adaptIntegrateBallSplit(f, n, xstar, width = 0, lowerLimit = rep(0, n - 1), 
    upperLimit = c(rep(pi, n - 2), 2 * pi), R = c(0, 1), tol = 1e-05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{Integrand function f(x)=f(x[1],...,x[n]).
}
  \item{n}{dimension of the space.  The sphere is an (n-1) dimensional
   manifold inside n-space, the ball is an n-dimensional solid.
}
  \item{lowerLimit}{Polar angular coordinates of lower limit
}
  \item{upperLimit}{Polar angular coordinates of upper limit
}
  \item{tol}{tolerance, the desired accuracy of the result.  The functions try to get abs(exact-value) < tol
}
  \item{\dots}{optional arguments passed to f.  If used, these should be specified with a tag, e. g. param1=7
}
  \item{R}{a numeric vector of length 2, integration is performed over the region with R[1] < radius < R[2].
}
  \item{xstar}{(n x m) matrix whose columns give the directions where the integrand changes quickly, where the 
  region will be subdivided to focus on that region.  (The length of a column vector is not used, just it's direction.)
}
  \item{width}{width of 'splitting regions', a vector of length m. If it is of length 1, 
  then that valued is repeated for each j in 1:m.
   If width[j]=0, the angular region is split just at the points given by the columns xstar[ ,j].
  If width[j] > 0, then angular region is split at an angle plus and minus width[j].
}
}
\details{Approximate the integral of f(x) over (part of) the sphere or ball in n-space.
The approach is simplistic: reparameterize the
region in polar coordinates.  For the sphere, this makes the region of integration a hyper-rectangle
in dimension (n-1) in the angle space (here the radius is fixed: R=1).  For the ball, the polar representation in terms
of angles and radius gives an region of integration that is an n dimensional hyper-rectangle.

The region of integration can be a subset of the sphere/ball by specifying a
patch/sector in polar coordinates.  
To integrate over a subregion, bounds for the polar integration have be specified.  For example,
in two dimensions, you can integrate over the top half of the circle by specifying lowerLimit=0.0
and upperLimit=pi to adaptIntegrateSphere.  Likewise for the ball, to integrate over the part of the annulus
 with inner radius .2 and outer radius .7 that is in the first quadrant, specify
lowerLimit=0.0, upperLimit=pi/2, R=c(.2,.7).
}
\value{For adaptIntegrateSphere and adaptIntegrateBall, the function returns a
 list containing several fields.  There is always a field
   \item{value}{Giving the approximate value of the integral.}
The other fields depend on the dimension: when n=2, the other fields are
what is returned by the function integrate( ) in base \R;
when n > 2, the other fields are the fields returned by function adaptIntegrate( ) in package cubature.

For adaptIntegrateSphereSplit and adaptIntegrateBallSplit, a single value is returned.  (This is 
because these functions make multiple calls to the adaptive integration routine and the 
results of each call are not saved.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{polar2rect}}, \code{\link{rect2polar}}
}
\examples{
f1 <- function( x ) { return(x[1]^2+3*x[2]+exp(x[3])) }
n <- 3
adaptIntegrateSphere( f1, n )
adaptIntegrateSphereSplit( f1, n, xstar=matrix(c(1,1,1),nrow=3) )
adaptIntegrateBall( f1, n )
adaptIntegrateBallSplit( f1, n, xstar=matrix(c(1,1,1),nrow=3) )

# test of adaptive integration with deliberate splitting
# function f3 has a sharp spike in the direction (1,2),
# elsewhere it has value 1
f3 <- function( x ) {
  x0 <- c(1.0,2.0)/sqrt(5.0)
  dist <- sqrt(sum( (x-x0)^2) )
  y <- 10-5000*dist
  y <- 1 + max(y,0)
  return(y)  }

# no splitting: this straightforward attempt at integration misses
# the spike and sees the integrand as =1 everywhere, so returns the arclengh 2*pi
n <- 2
adaptIntegrateSphere( f3, n ) 

# deliberate splitting at specified points, but still misses spike
# default width=0 splits the region of integration from [0,2*pi] to [0,a] U [a,2*pi],
# where tan(a)=2/1.
xstar <- matrix( c(1.0,2.0,-1.0,1.0), nrow=2 )
adaptIntegrateSphereSplit( f3, n, xstar=xstar )

# deliberate splitting around specified points, 'smart' choice of width gets the spike
# Here the region of integration is split into [0,a-.01] U [a-.01,a+.01] U [a+.01,2*pi]
adaptIntegrateSphereSplit( f3, n, xstar=xstar, width=c(0.01,0.01) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cubature }
\keyword{ multivariate integration }
