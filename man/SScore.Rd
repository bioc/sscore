\name{SScore }
\alias{SScore}
\title{ Compute S-Score values }
\description{
  Computes the S-Score values for a pair of Affymetrix GeneChips
}
\usage{
SScore(afbatch = stop("No CEL files specified"), classlabel = c(0,1), SF = NULL, 
       SDT = NULL, rm.outliers = TRUE, rm.mask = TRUE, rm.extra = NULL, 
       digits = NULL, verbose = FALSE, celfile.path = NULL, celfile.names = NULL)
}
\arguments{
  \item{afbatch}{ An \code{AffyBatch} object }
  \item{classlabel}{ A vector identifying the class for each column of the \code{AffyBatch} object }
  \item{SF}{ a list of Scale Factor (SF) values for each GeneChip }
  \item{SDT}{ a list of Standard Difference Threshold (SDT) values for each GeneChip }
  \item{rm.outliers}{ should the spots marked as 'OUTLIERS' be excluded from S-Score calculation? }
  \item{rm.mask}{ should the spots marked as 'MASKS' be excluded from S-Score calculation? }
  \item{rm.extra}{ if \code{TRUE}, overrides what is in \code{rm.mask} and \code{rm.outliers} }
  \item{digits}{ number of significant digits for S-Score values }
  \item{verbose}{ logical value. If \code{TRUE} it provides more detail of the S-Score calculations. }
  \item{celfile.path}{ character denoting the path for the *.CEL files corresponding to \code{afbatch} }
  \item{celfile.names}{ optional character vector containing the names of the *.CEL files }
}
\details{
  Computes S-Score values as described by Zhang et al. (2002).  \code{SScore} provides a simpler
  interface for comparing only two classes of GeneChips, while \code{SScoreBatch} compares 
  multiple pairs of chips.

  The classlabel consists of a vector with one entry for each column of the \code{AffyBatch} object.
  Each entry consists of a 0 or a 1 to identify the class to which the chip for the corresponding column
  belongs.  \code{SScore} will conduct a two-class test comparing all chips labeled 0 to all chips
  labeled 1.  If classlabel is not specified, it defaults to a two-chip comparison, compatible with
  previous versions of \code{SScore}.

  The SF and SDT factors are required for all calculations.  If \code{NULL}, these values will be
  calculated according to the Affymetrix Statistical Algorithms Description Document.  \code{digits}
  allows the specification of the number of significant digits for the S-Score values; if \code{NULL},
  the maximum number of significant digits are retained.
}
\value{
  An \code{ExpressionSet} with S-Score values in the \code{exprs} slot.
}
\references{ 

Zhang, L., Wang, L., Ravindranathan, A., Miles, M.F. 	(2002) A new algorithm for
 analysis of oligonucleotide arrays: application to expression profiling in mouse
 brain regions.  Journal of Molecular Biology, 317(2), pp. 225--35

Kerns, R.T., Zhang, L., Miles, M.F.  (2003) Application of the S-score algorithm for 
analysis of oligonucleotide microarrays.  Methods, 31(4), pp. 274--81
}
\author{ Richard Kennedy \email{rkennedy@vcu.edu}}
\note{ Based on C++ code by Li Zhang and Delphi code by Robnet Kerns }
\seealso{ \code{\link{SScoreBatch}},\code{\link{computeSFandSDT}},\code{\link{computeOutlier}} }
\examples{
  if (length(dir(pattern=".cel$")) != 0) {

  ## Read in the *.CEL files
  abatch <- ReadAffy()

  ## default calling method
  SScores <- SScore(abatch)

  ## specifying SF and SDT (gives same results as above)
  SfSdt <- computeSFandSDT(abatch)
  SScores <- SScore(abatch,SF=SfSdt$SF,SDT=SfSdt$SDT)

  ## specifying outlier and masked values should be included in calculations
  SScores <- SScore(abatch,rm.outliers=FALSE,rm.mask=FALSE)

  ## round results to 3 significant digits
  SScores <- SScore(abatch,digits=3)

  ## show verbose output
  SScores <- SScore(abatch,verbose=TRUE)

}

}
\keyword{ manip }
