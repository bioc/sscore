\name{SScoreBatch }
\alias{SScoreBatch}
\title{ Compute S-Score values }
\description{
  Computes the S-Score values for multiple pairs of Affymetrix GeneChips
}
\usage{
SScoreBatch(afbatch, compare, SF = NULL, SDT = NULL, 
            rm.outliers = TRUE, rm.mask = TRUE, rm.extra = TRUE, 
            digits = NULL, verbose = FALSE, celfile.path = NULL)
}
\arguments{
  \item{afbatch}{ An \code{AffyBatch} object }
  \item{compare}{ A matrix describing which chips to compare }
  \item{SF}{ a list of Scale Factor (SF) values for each GeneChip }
  \item{SDT}{ a list of Standard Difference Threshold (SDT) values for each GeneChip }
  \item{rm.outliers}{ should the spots marked as 'OUTLIERS' be excluded from S-Score calculation? }
  \item{rm.mask}{ should the spots marked as 'MASKS' be excluded from S-Score calculation? }
  \item{rm.extra}{ if \code{TRUE}, overrides what is in \code{rm.mask} and \code{rm.outliers} }
  \item{digits}{ number of significant digits for S-Score values }
  \item{verbose}{ logical value. If \code{TRUE} it provides more detail of the S-Score calculations. }
  \item{celfile.path}{ character denoting the path for the *.CEL files corresponding to \code{afbatch} }
}
\details{
  Computes S-Score values as described by Zhang et al. (2002).  \code{SScoreBatch} allows comparison
  of multiple pairs of chips, while \code{SScore} provides a simpler interface when comparing only
  two GeneChips.
  
  \code{compare} specifies how the pairwise comparisons are performed.  It is an N x 2 matrix,
  where N is the number of pairwise comparisons; each row of the matrix contains index in the \code{AffyBatch} object
  for the chips to be compared.  For example,

  \preformatted{
     1     3 
     4     2 
     5     9 
    10     2
     5     7
  }

  would do a comparison of chip 1 to chip 3, a comparison of chip 4 to chip 2, a comparison of chip 5
  to chip 9, and so on.  The columns in \code{exprSet} correspond to the rows of \code{compare}, so
  that the results of the first comparison are in column 1, the results of the second comparison are
  in column 2, and so on.

  The SF and SDT factors are required for all calculations.  If \code{NULL}, these values will be
  calculated according to the Affymetrix Statistical Algorithms Description Document.  
  probes.  \code{digits} allows the specification of the number of significant digits for the S-Score
  values; if \code{NULL}, the maximum number of significant digits are retained.
}
\value{
  An \code{exprSet} with S-Score values in the \code{exprs} slot and
  CorrDiff values in the \code{se.exprs}.
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
\seealso{ \code{\link{SScore}},\code{\link{computeSFandSDT}},\code{\link{computeOutlier}} }
\examples{
  if (length(dir(pattern=".cel$")) != 0) {

  ## Read in the *.CEL files
  abatch <- ReadAffy()

  ## default calling method
  SScores <- SScoreBatch(abatch)

  ## specifying SF and SDT (gives same results as above)
  SfSdt <- computeSFandSDT(abatch)

  ## specifying outlier and masked values should be included in calculations
  SScores <- SScoreBatch(abatch,SF=SfSdt$SF,SDT=SfSdt$SDT)
  SScores <- SScoreBatch(abatch,rm.outliers=FALSE,rm.mask=FALSE)

  ## round results to 3 significant digits
  SScores <- SScoreBatch(abatch,digits=3)

  ## show verbose output
  SScores <- SScoreBatch(abatch,verbose=TRUE)

}

}
\keyword{ manip }
