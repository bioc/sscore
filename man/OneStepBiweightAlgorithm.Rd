\name{OneStepBiweightAlgorithm }
\alias{OneStepBiweightAlgorithm}
\title{ One-step Tukey's biweight }
\description{
  Computes one-step Tukey's biweight on a vector.  Note that this implementation
  follows the Affymetrix code, which is different from the Tukey's biweight
  computed by the affy package.
}
\usage{
OneStepBiweightAlgorithm(x, c, epsilon)
}
\arguments{
  \item{x}{ vector of data }
  \item{c}{ tuning constant (see details) }
  \item{epsilon}{ fuzz value to avoid division by zero (see details) }
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
\keyword{ manip }