\name{SphericalMisc}
\alias{SphericalMisc}
\alias{adaptIntegrateSphereCheck}
\alias{nextGraySubset}
\alias{nextMultiIndex}
\alias{partitionRegion}
\alias{CheckUnitVectors}
\alias{Octants}
\alias{SubdivideSphereTriByOctant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Miscellaneous internal functions used by SphericalCubature package. 
}
\description{
These functions are not intended for general use, they are
only listed here to document their existence.

adaptIntegrateCheck is used by the adaptive integration functions
to check input parameters,
partitionRegion is used by the 'split' versions of the the adaptive
integration functions,  nextGraySubset is used by IntegrateSphereStroud11,
nextMultiIndex is used by adaptive integration functions.
}
\usage{
adaptIntegrateSphereCheck( n, lowerLimit, upperLimit, R, xstar, width )
partitionRegion( xstar, width, lowerLimit, upperLimit )
nextGraySubset( gray.list )
nextMultiIndex( j, size )
CheckUnitVectors( S, eps=1.0e-14 )
Octants( n, positive.only=FALSE )
SubdivideSphereTriByOctant( S, eps=1.0e-14 )
}
\arguments{
  \item{n}{dimension of the space}
  \item{lowerLimit}{lower angular limit for integration region}
  \item{upperLimit}{upper angular limit for integration region}
  \item{R}{inner and outer radii for integration region}
  \item{xstar}{directions where function changes rapidly}
  \item{width}{width of subdivisions}  
  \item{gray.list}{list used by Stroud integration}
  \item{j}{current multi-index}
  \item{size}{length of multi-index}
  \item{S}{a matrix or array specifying simplices}
  \item{eps}{tolerance used in test for unit vectors}
  \item{positive.only}{logical; if TRUE, only the first octant is returned, if FALSE, all 2^n octants are returned }
}

\details{
\code{Octant} returns an array of simplices in the V representation; each one gives the points on the axes that are on the boundary of one octant.
\code{CheckSphereTri} performs validity checks on the simplices in S.
}
