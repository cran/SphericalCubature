\name{integrateSphereStroud11}
\alias{integrateSphereStroud11}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Integrate a function over the sphere in n-dimensions.
}
\description{Approximate the integral of a function f(x)=f(x[1],...,x[n])
over the unit sphere in n-space using Stroud's method of degree 11.
}
\usage{
integrateSphereStroud11(f, n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{function f(x)=f(x[1],...,x[n]) to integrate
}
  \item{n}{dimension of the space, implemented for n in the range 3:16.
}
  \item{\dots}{optional arguments passed to f( ).  If these are specified, they should be labeled with a tag, 
  e.g. param1=3.4
}
}
\details{This method works if the integrand is smooth.
If the function changes rapidly, adaptive integration can be tried as described in 'See Also' below.
}
\value{A single number, the approximation to the integral.
}
\references{
Stroud integration and related functions, adapted from fortran code by John Burkhart found at \cr
  http://people.sc.fsu.edu/~jburkardt/f77_src/stroud/stroud.html \cr
Based on the book by A. H. Stroud, Approximate Calculation of
multiple integrals, 1971, page 296-297.
}
\seealso{
\code{\link{adaptIntegrateSpherePolar}}, \code{\link{adaptIntegrateBallPolar}}, \code{\link{adaptIntegrateSphereTri}}
}
\examples{
f2 <- function( x ) { return(x[1]^2) }
integrateSphereStroud11( f2, n=3 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate integration }
\keyword{ cubature }
