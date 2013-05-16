\name{SphericalMisc}
\alias{SphericalMisc}
\alias{adaptIntegrateCheck}
\alias{nextGraySubset}
\alias{nextMultiIndex}
\alias{partitionRegion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Miscellaneous functions used by SphericalCubature package. 
}
\description{
These functions are not intended for general use, they are
only documented here to acknowledge their existence and to get
\R build system from complaining about undocumented functions.

adaptIntegrateCheck is used by the adaptive integration functions
to check input parameters,
partitionRegion is used by the 'split' versions of the the adaptive
integration functions,  nextGraySubset is used by IntegrateSphereStroud11,
nextMultiIndex is used by adaptive integration functions.
}
\usage{
adaptIntegrateCheck( n, lowerLimit, upperLimit, R, xstar, width )
partitionRegion( xstar, width, lowerLimit, upperLimit )
nextGraySubset( gray.list )
nextMultiIndex( j, size )
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
}
