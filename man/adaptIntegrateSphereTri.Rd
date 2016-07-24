\name{adaptIntegrateSphereTri}
\alias{adaptIntegrateSphereTri}
\alias{adaptIntegrateSphereTri3d}
\alias{adaptIntegrateSphereTriI0}
\alias{adaptIntegrateSphereTriI1}
\alias{adaptIntegrateSphereTriSubdivideK}
\title{Adaptive integration over spherical triangles}
\description{Adaptively integrate a function over a set of spherical triangles.
Function \code{adaptIntegrateSphereTri} uses hyperspherical triangles and 
works in n-dimensions; it uses function \code{adaptIntegrateSimplex} in package SimplicialCubature, which is 
based on code of Alan Genz.  \code{adaptIntegrateSphereTri3d} works only in 3-dimensions and is 
described in the paper by N. Boal and F-J. Sayas; it is not as sophisticated and is slower than \code{adaptIntegrateSphereTri}, but can be useful and is self contained.  
}
\usage{
adaptIntegrateSphereTri( f, S, fDim=1L, maxEvals=20000L, absError=0.0, 
    tol=1.0e-5, integRule=3L, partitionInfo=FALSE, ...  ) 
adaptIntegrateSphereTri3d( f, S, maxRef=50, relTol=0.001, maxTri=50000, gamma=1.5 ) 
adaptIntegrateSphereTriI0( f, K )
adaptIntegrateSphereTriI1( f, K )
adaptIntegrateSphereTriSubdivideK( K ) 
}
\arguments{
  \item{f}{function f defined on the sphere}
  \item{S}{array of hyperspherical triangles, dim(S)=c(n,n,nS).  Columns of S should be points on the unit sphere: sum(S[,i,j]^2)=1. It is required that every simplex S[,i,j] is 
  contained within a single octant; it can be any octant, but not simplices can't have vertices in different octants.  This will happend automatically if Octants( ) is 
  used to generate octants, or if S is a tessellationg coming from function \code{UnitSphere} in package \code{mvmesh} (or a subdivision of such a partition).}
  \item{fDim}{integer dimension of the integrand function}
  \item{maxEvals}{maximum number of evaluations allowed}
  \item{absError}{desired absolute error}
  \item{tol}{desired relative tolerance}
  \item{integRule}{integration rule to use in call to function \code{adsimp}}
  \item{partitionInfo}{if TRUE, return the final partition after subdivision}
  \item{...}{optional arguments to function f(x,...)}  
  \item{maxRef}{maximum number of refinements allowed}
  \item{relTol}{desired relative tolerance}
  \item{maxTri}{maximum number of triangles allowed while refining}
  \item{gamma}{threshold parameter for when a triangle is subdivided, should be >= 1}
  \item{K}{a single spherical triangle in 3 space, dim(K)=c(3,3)}
}
\details{\code{adaptIntegrateSphereTri3d} takes as input a function 
f defined on the unit sphere in 3-dimensions and a list of spherical triangles K0 and attempts to 
integrate f over the part of the unit sphere described by K0.  The triangles in K0 should individually contained 
in a hemisphere. The algorithm estimates the integral over each triangle, then estimates error on each triangle, 
and subdivides triangles with large errors.  This is repeated until either the desired accuracy is achieved, 
or there are too many subdivisions.

Functions \code{adaptIntegrateSphereI0}, \code{adaptIntegrateSphereI1}, and \code{adaptIntegrateSphereSubdivideK} are internal 
functions that compute respectively I0 for a spherical triangle, I1 for a spherical triangle, and subdivide
spherical triangle K into 4 smaller spherical triangles (using midpoints of each side).

\code{adaptIntegrateSphereTri}  takes as input a function 
f defined on the unit sphere in n-dimensions and a list of hyperspherical triangles S and attempts to 
integrate f over (part of) the unit sphere described by S.  It uses the \R package \code{SimplicialCubature} to
evaluate the integrals.  The hyperspherical triangles in S should individually be contained 
in a hemisphere. The algorithm estimates the integral over each triangle, then estimates error on each triangle, 
and subdivides triangles with large errors.  This is repeated until either the desired accuracy is achieved, 
or there are too many subdivisions.  This function is more general than \code{adaptIntegrateSphereTri3d} in two 
ways: (1) it works in dimension n >= 2, and (2) it allows vector integrands f.  It also is generaly faster
than \code{adaptIntegrateSphereTri}. Use the function \code{Octants} to get a rough triangulation of
the sphere; see the examples below.  For finer resolution triangulation can be obtained from 
\code{UnitSphere} in \R package \code{mvmesh}.    
For description of \code{adaptIntegrateSphereTri3d} see www.unizar.es/galdeano/actas_pau/PDFVIII/pp61-69.pdf. 


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
# test of polynomial integration f(s)=s[1]^2
f <- function( s ) { return( s[1]^2 ) }
n <- 3

# integrate over whole sphere
S <- Octants( n )
a <- adaptIntegrateSphereTri( f, S )
b <- adaptIntegrateSphereTri3d( f, S )
# exact answer, adaptIntegrateSphereTri approximation, adaptIntegrateSphereTri3d approximation
sphereArea(n)/n; a$integral; b$integral

# integrate over first octant only
S <- Octants( n, positive.only=TRUE ) 
a <- adaptIntegrateSphereTri( f, S )
b <- adaptIntegrateSphereTri3d( f, S )
# exact answer, adaptIntegrateSphereTri approximation, adaptIntegrateSphereTri3d approximation
sphereArea(n)/(8*n); a$integral; b$integral


}


