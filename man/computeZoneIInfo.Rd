\name{computeZoneIInfo }
\alias{computeZoneIInfo}
\title{ Compute zone background and noise }
\description{
  Computes the background and noise for a given zone
  of a single Affymetrix GeneChip
}
\usage{
computeZoneIInfo(ZoneInfo, NumberBGCells)
}
\arguments{
  \item{ZoneInfo}{ vector of intensities in a given zone }
  \item{NumberBGCells}{ number of background cells for the GeneChip }
}
\details{
  Calculates background and noise for a zone using the algorithms described in the
  Affymetrix Statistical Algorithms Description Document.  This is an internal
  function that will generally not be accessed directly.
}
\value{\code{computeZoneIInfo} returns a list containing the following components:
  \item{background }{ background value for the given zone }
  \item{noise }{ noise value for the given zone }
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
