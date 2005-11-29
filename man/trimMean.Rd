\name{trimMean}
\alias{trimMean}
\title{ Compute trimmed mean for a vector }
\description{
  Computes the trimmed mean for a vector.  Note that this implementation
  follows the Affymetrix code, which gives different results than the
  standard R function mean().
}
\usage{
trimMean(vec, p1, p2)
}
\arguments{
  \item{vec}{ vector of values }
  \item{p1}{ lower percentage for trimming }
  \item{p2}{ upper percentage for trimming }
}
\details{
  The details can be found in the given reference.
}
\value{
  A numeric value
}
\references{ 

Affymetrix (2002) Statistical Algorithms Description Document, Affymetrix
Inc., Santa Clara, CA,
whitepaper.
\url{http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf}
}
\author{ Richard Kennedy \email{rkennedy@vcu.edu}}
\note{ Based on Affymetrix MAS5 Statistical SDK source code 
\url{http://www.affymetrix.com/Auth/support/developer/stat_sdk/STAT_SDK_source.zip}
}
\keyword{ univar }