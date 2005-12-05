\name{computeSFandSDT }
\alias{computeSFandSDT}
\title{ Compute SF and SDT values }
\description{
  Computes the scaling factor (SF) and statistical difference threshold (SDT) values
  of Affymetrix GeneChips, for use in calculating S-Score values
}
\usage{
computeSFandSDT(afbatch, TGT = 500, digits = NULL, verbose = FALSE, plot.histogram = FALSE, celfile.path = NULL)
}
\arguments{
  \item{afbatch}{ An \code{AffyBatch} object }
  \item{TGT}{ the target intensity to which the arrays should be scaled }
  \item{digits}{ number of significant digits for SF and SDT values }
  \item{verbose}{ logical value. If \code{TRUE} it provides more detail of the SF and SDT calculations. }
  \item{plot.histogram}{ logical value.  if \code{TRUE} it plots a histogram of intensities }
  \item{celfile.path}{ character denoting the path for the *.CEL files corresponding to \code{afbatch} }
}
\details{
  Calculates SF and SDT factors using the algorithms described in the Affymetrix Statistical
  Algorithms Description Document.  The SF and SDT may be used in the calculation of S-Score
  values, or may be useful in their own right.  One SF and SDT value is calculated for each
  GeneChip, which are arranged in the same order as the columns in the \code{AffyBatch} object.
}
\value{\code{computeSFandSDT} returns a list containing the following components:
  \item{SF }{ SF values, one for each GeneChip }
  \item{SDT }{ SDT values, one for each GeneChip }
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
\examples{
  if (length(dir(pattern=".cel$")) != 0) {

  ## Read in the *.CEL files
  abatch <- ReadAffy()

  ## compute SF and SDT
  SfSdt <- computeSFandSDT(abatch)

  ## show verbose output
  SfSdt <- computeSFandSDT(abatch,verbose=TRUE)

  ## plot PM and MM histograms for each *.CEL file
  SfSdt <- computeSFandSDT(abatch,plot.histogram=TRUE)

}

}
\keyword{ manip }
