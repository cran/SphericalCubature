\name{adaptIntegrateBallTri}
\alias{adaptIntegrateBallTri}
\alias{adaptIntegrateBallRadial}
\title{Adaptive integration over the unit ball}
\description{Adaptively integrate a function over the ball, specified by a set of spherical triangles
to define a ball (or a part of a ball).
Function \code{adaptIntegrateBallTri()} uses spherical triangles and 
works in n-dimensions; it uses function \code{adaptIntegrateSimplex()} in \R package \code{SimplicialCubature}, which is 
based on code of Alan Genz.  \code{adaptIntegrateBallRadial()} integrates radial functions of the form 
f(x) = g(|x|) over the unit ball.
}
\usage{
adaptIntegrateBallTri( f, n, S=Orthants(n), fDim=1L, maxEvals=20000L, absError=0.0, 
    tol=1.0e-5, integRule=3L, partitionInfo=FALSE, ...  ) 
adaptIntegrateBallRadial( g, n, fDim=1, maxEvals=20000L, absError=0.0, 
    tol=1e-05, integRule=3L, partitionInfo=FALSE, ... )
}
\arguments{
  \item{f}{integrand function f defined on the sphere in R^n}
  \item{g}{inegrand function g defined on the real line}
  \item{n}{dimension of the space}
  \item{S}{array of spherical triangles, dim(S)=c(n,n,nS).  Columns of S should be points on the unit sphere: sum(S[,i,j]^2)=1. Execution will be faster if every simplex S[,,j] is 
  contained within any single orthant.  This will happend automatically if function \code{Orthants}  is 
  used to generate orthants, or if S is a tessellation coming from function \code{UnitSphere} in package \code{mvmesh}. If one or more simplices
  interect multiple orthants, the simplices will automatically be subdivided so that each subsimplex is in a single orthant.}
  \item{fDim}{integer dimension of the integrand function}
  \item{maxEvals}{maximum number of evaluations allowed}
  \item{absError}{desired absolute error}
  \item{tol}{desired relative tolerance}
  \item{integRule}{integration rule to use in call to function \code{adsimp}}
  \item{partitionInfo}{if TRUE, return the final partition after subdivision}
  \item{...}{optional arguments to function f(x,...) or g(x,...)}  
}
\details{\code{adaptIntegrateBallTri()}  takes as input a function 
f defined on the unit sphere in n-dimensions and a list of spherical triangles S and attempts to 
integrate f over (part of) the unit sphere described by S.  It uses the \R package \code{SimplicialCubature} to
evaluate the integrals.  The spherical triangles in S should individually be contained 
in an orthant.   

If the integrand is nonsmooth, you can specify a set of spherical triangles that focus the cubature routines 
on that region.  See the example below with integrand function \code{f3}.
}

\value{A list containing 
\item{status}{a string describing result, ideally it should be "success", otherwise it is an error/warning message.}
\item{integral}{approximation to the value of the integral}
\item{I0}{vector of approximate integral over each triangle in K}
\item{numRef}{number of refinements}
\item{nk}{number of triangles in K}
\item{K}{array of spherical triangles after subdivision, dim(K)=c(3,3,nk)}
\item{est.error}{estimated error}
\item{subsimplices}{
      if partitionInfo=TRUE, this gives an array of subsimplices, see function adsimp for more details.}
\item{subsimplicesIntegral}{
      if partitionInfo=TRUE, this array gives estimated values of each component of the integral on each 
            subsimplex, see function adsimp for more details.}
\item{subsimplicesAbsError}{
      if partitionInfo=TRUE, this array gives estimated values of the absolute error of each component of the integral on each 
            subsimplex, see function adsimp for more details.}
\item{subsimplicesVolume}{
      if partitionInfo=TRUE, vector of m-dim. volumes of subsimplices; this is not d-dim. volume if m < n.}

}
\examples{
# integrate over ball in R^3
n <- 3
f <- function( x ) { x[1]^2  }
adaptIntegrateBallTri( f, n )

# integrate over first orthant only
S <- Orthants( n, positive.only=TRUE ) 
a <- adaptIntegrateSphereTri( f, n, S )
b <- adaptIntegrateSphereTri3d( f, S )
# exact answer, adaptIntegrateSphereTri approximation, adaptIntegrateSphereTri3d approximation
sphereArea(n)/(8*n); a$integral; b$integral

# integrate a vector valued function
f2 <- function( x ) { c(x[1]^2,x[2]^3) }
adaptIntegrateBallTri( f2, n=2, fDim=2 )

# example of specifiying spherical triangles that make the integration easier
f3 <- function( x ) { sqrt( abs( x[1]-3*x[2] ) ) } # has a cusp along the line x[1]=3*x[2]
a <- adaptIntegrateBallTri( f3, n=2, absError=0.0001, partitionInfo=TRUE )
str(a) # note that returnCode = 1, e.g. maxEvals exceeded

# define problem specific spherical triangles, using direction of the cusp
phi <- atan(1/3); c1 <- cos(phi); s1 <- sin(phi)
S <- array( c(1,0,c1,s1,  c1,s1,0,1,    0,1,-1,0,   -1,0,-c1,-s1,  
            -c1,-s1,0,-1, 0,-1,1,0), dim=c(2,2,6) )
b <- adaptIntegrateBallTri( f3, n=2, S, absError=0.0001, partitionInfo=TRUE )
str(b) # here returnCode=0, less than 1/2 the function evaluations and smaller estAbsError

# integrate x[1]^2 over nested balls of radius 2 and 5 (see discussion in ?SphericalCubature)
f4 <- function( x ) { c( 2^2 * (2*x[1])^2, 5^2 * (5*x[1])^2 ) }
bb <- adaptIntegrateBallTri( f4, n=2, fDim=2 )
str(bb)
bb$integral[2]-bb$integral[1] # = integral of x[1]^2 over the annulus 2 <= |x| <= 5

}
