\name{sphereArea}
\alias{sphereArea}
\alias{ballVolume}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Surface area of spheres, volumes of balls in n-dimensions.
}
\description{Calculates the (n-1) dimensional surface area 
of a sphere and the n dimensional volume of a ball in n-space.
}
\usage{
sphereArea(n, R = 1)
ballVolume(n, R = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Dimension of the space.
}
  \item{R}{Radius of the sphere/ball.}
}
\value{Single number that is the area of the sphere/volume of the ball.
}

\examples{
sphereArea(n=5)
ballVolume(n=5)
}
