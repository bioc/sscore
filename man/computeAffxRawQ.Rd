\name{computeAffxRawQ }
\alias{computeAffxRawQ}
\title{ Computes RawQ using affxparser routines}
\description{
  Computes the RawQ value of a single Affymetrix GeneChip
}
\usage{
computeAffxRawQ(intensity, stdvs, pixels, probe.index, probe.zoneID, bgCells, NumberZones)
}
\arguments{
  \item{intensity}{ vector of intensities for the GeneChip }
  \item{stdvs} { matrix of standard deviations for probe intensities from *.CEL file }
  \item{pixels} { matrix of number of pixels for probe intensities from *.CEL file }
  \item{probe.index}{ vector of indices for each probe }
  \item{probe.zoneID}{ vector of zone ID numbers for each probe }
  \item{bgCells}{ number of background cells for the GeneChip }
  \item{NumberZones}{ number of zones on the GeneChip }
}
\details{
  Calculates RawQ using the algorithms described in the Affymetrix Statistical Algorithms
  Description Document for a single GeneChip.  This is an internal function that will
  generally not be accessed directly.
}
\value{
  the RawQ value for the given array
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
