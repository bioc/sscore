\name{computeOutlier }
\alias{computeOutlier}
\title{ Compute outlier probes  }
\description{
  Computes outlier and/or mask probes for a set of Affymetrix GeneChips
  that will be excluded from S-Score calculations
}
\usage{
computeOutlier(afbatch, rm.mask = TRUE, rm.outliers = TRUE, rm.extra = TRUE, celfile.path = NULL)
}
\arguments{
  \item{afbatch}{ An \code{AffyBatch} object }
  \item{rm.mask}{ should the spots marked as 'MASKS' be excluded from S-Score calculation? }
  \item{rm.outliers}{ should the spots marked as 'OUTLIERS' be excluded from S-Score calculation? }
  \item{rm.extra}{ if \code{TRUE}, overrides what is in \code{rm.mask} and \code{rm.outliers} }
  \item{celfile.path}{ character denoting the path for the *.CEL files corresponding to \code{afbatch} }
}
\details{
  Computes the outlier and / or mask probes for an \code{AffyBatch} object.  These are
  returned in matrix form, with one probe per row and one chip per column.  The value
  of each location in the matrix will be \code{TRUE} if the corresponding probe is an outlier /
  masked value and \code{FALSE} if it is not.  The options may be set to exclude only outlier values,
  only mask values, or both.  The probes are be arranged in the same row order as
  the intensity values.  Note that this function assumes the *.CEL files are still available
  in the directory given by \code{celfile.path} (or the current directory if \code{celfile.path} is not specified).
}
\value{
  a matrix containing the list of outliers / masked values for the given \code{AffyBatch} object.
}
\author{ Richard Kennedy \email{rkennedy@vcu.edu}}
\note{ Based on C++ code by Li Zhang and Delphi code by Robnet Kerns }
\examples{
  if (length(dir(pattern=".cel$")) != 0) {
  abatch <- ReadAffy()
  outlier <- computeOutlier(abatch)
}
}
\keyword{ manip }
