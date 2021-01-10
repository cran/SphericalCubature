\name{integrateSpherePolynomial}
\alias{integrateSpherePolynomial}
\alias{integrateBallPolynomial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Integration of polynomials over sphere or ball.}
\description{Exact integration of polynomial over sphere or ball in n-dimensions.
}
\usage{
integrateSpherePolynomial(p, valueOnly = TRUE)
integrateBallPolynomial(p, R = c(0, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a list specifying the coefficients and powers of the polynomial.  See details below
}
  \item{valueOnly}{boolean saying whether to return only the value of the integral, or return both
  the value and a intermediate terms.  These intermediate terms are used by integrateBallPolynomial( ).
}
  \item{R}{inner and outer radius of the annular region: R[1] <= radius <= R[2].
}
}
\details{
Compute the exact integral over the sphere in n dimensions of
a polynomial \cr
        p(x[1],...,x[n])=sum (coef[i] * x[1]^k[i,1] * ... * x[n]^k[i,n]), \cr
where the sum is over i=1,...,m.  The polynomial is specified as a list \code{p}
with fields
\itemize{
   \item{coef, an m-vector of doubles}
   \item{k, an (m x n) matrix of integers}}
   m and n are given implicitly in the sizes of these arrays
output is normally just a number, the value of the integral.  If  
integrateSpherePolynomial is called with valueOnly=FALSE, a list with two fields:
\itemize{  \item{integral, a double containing the value of the integral}
   \item{term, a vector of length m of values used in function IntegratePolynomialBall( )}}
}
\value{\code{integrateSpherePolynomial()}  normally just returns a value of the integral, but if 
 valueOnly=FALSE, it will return a list containing the value and intermediate terms.  These
 intermediate terms correspond to the integral of each monomial term in the polynomial; they
 are used by \code{integrateBallPolynomial()}.  
 
 integrateBallPolynomial( ) returns just the value of the integral.

}
\references{Method is from How to Integrate a Polynomial over a Sphere, 
by G. Folland (2001), MAA Monthly 108, pg. 446-448.
}

\examples{
n <- 3
# specify the polynomial p(x) = 1.0 * x[1]^2 * x[2]^0 * x[3]^0 + 7.0 * x[1]^0 * x[2]^3 * x[3]
p <- list(coef=c(1.0,7.0),k=matrix( c(2L,0L,0L,0L,3L,0L), byrow=TRUE, nrow=2) )
integrateSpherePolynomial( p )
integrateBallPolynomial( p )

# compare to adaptive integration
f4 <- function( x ) { return( x[1]^2 + 7*x[2]^2*x[3] ) }
adaptIntegrateSpherePolar( f4, n )$integral
adaptIntegrateBallPolar( f4, n )$integral
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cubature }
\keyword{ multivariate integration }
