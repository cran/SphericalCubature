\name{rect2polar}
\alias{rect2polar}
\alias{polar2rect}
\title{n-dimensional polar coordinate transformations}
\description{Convert between polar and rectangular coordinates in n-dimensions.
The point  (x[1],...,x[n]) in rectangular coordinates corresponds to the point  (r,phi[1],...,phi[n-1])
in polar coordinates.  
}
\usage{
polar2rect(r, phi)
rect2polar(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{a vector of radii of length m.
}
  \item{phi}{a (n-1) x m matrix of angles.
}
  \item{x}{ (n x m) matrix, with column j being the point in n-dimensional space.
}
}
\details{n dimensional polar coordinates are given by the following: \cr
 rectangular x=(x[1],...,x[n])  corresponds to
   polar  (r,phi[1],...,phi[n-1]) by \cr
 x[1]  = r*cos(phi[1]) \cr
 x[2]  = r*sin(phi[1])*cos(phi[2]) \cr
 x[3]  = r*sin(phi[1])*sin(phi[2])*cos \cr
   ... \cr
 x[n-1]= r*sin(phi[1])*sin(phi[2])*...*sin(phi[n-2])*cos(phi[n-1]) \cr
 x[n]  = r*sin(phi[1])*sin(phi[2])*...*sin(phi[n-2])*sin(phi[n-1]) \cr
 
 Here phi[1],...,phi[n-2] in [0,pi), and phi[n-1] in [0,2*pi).
 For multivariate integration, the Jacobian of the above tranformation
  is J(phi) = r^(n-1) * prod( sin(phi[1:(n-2)])^((n-2):1) );
  note that phi[n-1] does not appear in the Jacobian.
}
\value{For polar2rect( ), an (n x m) matrix of rectangular coordinates.

 For rect2polar( ), a list with fields:
 \item{r}{a vector of length m containing the radii}
 \item{phi}{an (n x m) matrix of angles}
}


\examples{
x <- matrix( 1:9, nrow=3 )
x
a <- rect2polar( x )
a
polar2rect( a$r, a$phi )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ polar coordinates }
