\name{computeZonenoise }
\alias{computeZonenoise}
\title{ Computes zone noise }
\description{
  Computes the noise (average standard error) of the probe intensities
  for a single Affymetrix GeneChip
}
\usage{
computeZonenoise(index, intensity, stdv, npixels, bgCells)
}
\arguments{
  \item{index}{ vector of indices for probes in the given zone }
  \item{intensity}{ vector of intensities for the GeneChip }
  \item{stdv}{ vector of standard deviations for the GeneChip }
  \item{npixels}{ vector containing number of pixels for each probe of the GeneChip }
  \item{bgCells}{ number of background cells on the GeneChip }
}
\details{
  Calculates the noise (average standard error) of the probes in a given zone, using
  the algorithms described in the Affymetrix Statistical Algorithms Description Document,
  for a single GeneChip.  This is an internal function that will generally not be
  accessed directly.
}
\value{
  the noise of the probes for the given array
}
\references{ 

Affymetrix (2002) Statistical Algorithms Description Document, Affymetrix
Inc., Santa Clara, CA,
whitepaper.
\url{http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf}
}
\author{ Richard Kennedy \email{rkennedy@vcu.edu}}
\note{ Based on Affymetrix MAS5 Statistical SDK source code 
\url{http://www.affymetrix.com/Auth/support/developer/stat_sdk/STAT_SDK_source.zip},
as well as C++ code by Li Zhang and Delphi code by Robnet Kerns }
\keyword{ manip }
