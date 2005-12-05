SScoreBatch <- function(afbatch = stop("No CEL files specified"),
 compare =  stop("No list of comparisons given"),  SF = NULL,SDT =
 NULL, rm.outliers = TRUE,rm.mask = TRUE, rm.extra = TRUE, digits = 
 NULL,verbose = FALSE,celfile.path= NULL) {
#######################################################################
#
# This function computes the S-Score values for multiple probe sets in
# batch fashion.  It does this by doing consecutive pairwise
# comparisons specified by the user.
#
# Arguments:
#	afbatch - an AffyBatch object containing the data on the chips
#                 to be compared
#	compare - an N x 2 matrix, where N is the number of pairwise 
#                 comparisons; each row of the matrix contains the
#                 afbatch indices of the chips to be compared.  For
#                 example,
#				1	3
#				4	2
#				5	9
#				...
#
#		  would do a comparison of chip 1 to chip 3, a compar-
#                 ison of chip 4 to chip 2, a comparison of chip 5 to
#                 chip 9, etc.
#	SF     -  a vector of SF values (the Scale Factor, which is
#                 part of the general output from the Affymetrix GCOS
#                 software), one for each chip.  If SF = NULL, these
#                 will be calculated internally
#       SDT    -  a vector of SDT values (the Standard Difference
#                 Threshold, which can be derived from the general 
#                 output from the Affymetrix GCOS software using the
#                 formula SDT = 4 * RawQ * Scale Factor), one for each
#                 chip.  If SDT = NULL, these will be calculated
#                 internally
#	rm.outliers,
#       rm.mask,
#       rm.extra - logical value to remove outliers and/or masked 
#                  values.  If rm.extra = TRUE, it overrides the values
#                  of rm.outliers and rm.mask.  Works identically to
#                  these parameters in ReadAffy(), which it calls
#	digits - the number of significant digits to retain when
#                returning S-Score values
#	verbose  - logical value indicating whether additional informa-
#                  tion on the calculations should be displayed
#       celfile.path - character string giving the directory path to
#                      the *.CEL files being read, or NULL if the *.CEL
#                      files are in the current directory
#
# Value:
#	an object of class ExprSet, with the exprs slot containing the
#       S-score values and the se.exprs slot containing the CorrDiff
#       values
#
#######################################################################

# check the comparison matrix
	if (any(compare < 1) | any(compare > length(afbatch)))
		stop("Comparison index does not match number of CEL files")
	if (NCOL(compare) < 2)
		stop("Must have two chips for comparisons")
	if (ncol(compare) > 2)
		warning("More than two chips listed for each comparison.  Only the first two of each will be used")

# calculate SF and SDT, if not specified by the user, as these must
# always be available
	if (is.null(SF) | is.null(SDT))
		sfsdtlist <- computeSFandSDT(afbatch,digits=3,celfile.path=celfile.path)

	if (is.null(SF))
		SF <- sfsdtlist$SF
	if (any(SF <= 0)) 
		stop("SF values must be positive")
	if (length(SF) != length(afbatch))
		stop("Must be one SF value for each CEL file")

	if (is.null(SDT))
		SDT <- sfsdtlist$SDT
	if (any(SDT <= 0)) 
		stop("SDT values must be positive")
	if (length(SDT) != length(afbatch))
		stop("Must be one SDT value for each CEL file")

# calculate the outlier matrix, if excluding outliers / masked values
	if (any(rm.outliers,rm.mask,rm.extra))
		outlier <- computeOutlier(afbatch,rm.outliers=rm.outliers,rm.mask=rm.mask,rm.extra=rm.extra,celfile.path=celfile.path)

#######################################################################
# initialize variables needed in later calculations.  For the PM/MM
# index, these are calculated outside of the first loop, since they are
# constant for a given CDF.
	m.gamma <- 0.1
	probenames <- geneNames(afbatch)
	pmidx <- pmindex(afbatch)
	mmidx <- mmindex(afbatch)
	Score <- CorrDiff <- matrix(0.0,nrow=length(probenames),ncol=nrow(compare))

# this loops through each row of the compare matrix, which contains the
# indices of the pairwise comparisons; within the loop is the
# calculation for one pair of chips
	for (j in 1:nrow(compare)) {

# initialize variables needed in later calculations.  For idx1/idx2,
# these are assigned to reduce unwieldy expressions in subsequent code.
# For the intens1/intens2 and max1/max2 values, these are calculated
# outside of the second loop, since they are constant for a given chip.
		idx1 <- compare[j,1]
		idx2 <- compare[j,2]
		writeLines("Computing S-score values")
		intens1 <- intensity(afbatch[,idx1])*SF[idx1]
		intens2 <- intensity(afbatch[,idx2])*SF[idx2]
		max1 <- max(pm(afbatch[,idx1]),mm(afbatch[,idx1]))*SF[idx1]
		max2 <- max(pm(afbatch[,idx2]),mm(afbatch[,idx2]))*SF[idx2]

# this loops through each of the probesets on a given pair of chips
		for (i in 1:length(probenames)) {

# get the PM and MM values, as well as minimum intensity values, for
# the given probeset on each chip of the pair
			PM1 <- intens1[pmidx[[i]]]
			MM1 <- intens1[mmidx[[i]]]
			PM2 <- intens2[pmidx[[i]]]
			MM2 <- intens2[mmidx[[i]]]
			min1 <- min(PM1,MM1)
			min2 <- min(PM2,MM2)
	 
# adjust each of the PM and MM intensities relative to the minimum
# values
			PM1 <- PM1 - min1 
			PM2 <- PM2 - min2 
			MM1 <- MM1 - min1 
			MM2 <- MM2 - min2

# find the index of the probe pairs of the probeset to use in
# calculations.  A probe pair is used if it is not "saturated" (i.e.,
# the intensity is less than the maximum - minimum) and if it is not
# identified as an outlier / masked value in the .CEL file
			index <- (PM1<max1-min1)&(PM2<max2-min2)&(MM1<max1-min1)&(MM2<max2-min2)
			if (any(rm.outliers,rm.mask,rm.extra)) {
				outlier1 <- outlier[pmidx[[i]],idx1] | outlier[mmidx[[i]],idx1]
				outlier2 <- outlier[pmidx[[i]],idx2] | outlier[mmidx[[i]],idx2] 
				index <- index & (!outlier1) & (!outlier2)
			}
			N <- sum(as.integer(index))

# find the PM-MM differences for the probeset on each of the two chips in the pair
			diff1 <- (PM1-MM1)[index] 
			diff2 <- (PM2-MM2)[index] 

# this is the "raw" S-Score value (without N or alpha in calculation) -
# compare to S-Score calculation in J Mol Biol paper
			f.err <- (diff1-diff2)/sqrt(m.gamma*m.gamma*(diff1*diff1+diff2*diff2)+(SDT[idx1]*SDT[idx1]+SDT[idx2]*SDT[idx2])) 

# threshold or impute outlying f.err values.  The cutoff of 15 was
# arbitrarily decided in the original version; how it was determined
# is unknown
			f.err[f.err > 15.0] <- 15.0 
			f.err[f.err < -15.0] <- -15.0 
	 
			Score[i,j] <- sum(f.err) 
 
# estimate the variance / covariance values, for calculating the
# CorrDiff
			Sxx <- sum(diff1*diff1) 
			Syy <- sum(diff2*diff2) 
			Sxy <- sum(diff1*diff2) 
 
			Sx <- 0 
			Sy <- 0 

# transform the S-Score estimate by dividing by a function of the
# number of probes in the probeset
			if (N > 0) 
				Score[i,j] <- Score[i,j] / sqrt(N) else Score[i,j] <- 0 
 
# calculate the CorrDiff.  CorrDiffs below the threshold of 1e-3 are
# imputed to be 0, which was also arbitrarily decided in the first
# version
			if (N>2 && ((Sxx-Sx*Sx/N)*(Syy-Sy*Sy/N) > 1.e-3)) 
				CorrDiff[i,j] <- (Sxy-Sx*Sy/N)/sqrt((Sxx-Sx*Sx/N)*(Syy-Sy*Sy/N)) else CorrDiff[i,j] <- 0.0 
		}

# now renormalize the S-Score values, which gives alpha
		writeLines("Renormalizing S-scores")
		x <- Score[,j]
		Sx <- sum(x)
		Sxx <- sum(x*x)

# calculate the mean and standard deviation of the entire set of
# S-Scores
		Sstdev <- sqrt((Sxx-Sx*Sx/length(Score[,j]))/length(Score[,j]))
		meanSx <- Sx/length(Score[,j])

# find the trimmed S-score values, using a cutoff of those S-Scores
# within 3 standard deviations of the mean
		x <- Score[,j]-meanSx;
		x <- x[abs(x) < 3*Sstdev]
		Sx <- sum(x)
		Sxx <- sum(x*x)
		num <- length(x)

# calculate the trimmed mean and standard deviation.  Again, the cutoff
# of 0.01 was arbitrarily decided in the first version
		Sstdev <- ((Sxx-Sx*Sx/num)/num)
		if (Sstdev < 0.01) Sstdev <- 1.0 else Sstdev <- sqrt(Sstdev)
		m.alpha <- Sstdev
		meanSx <- Sx/num+meanSx

# perform the renormalization, using the trimmed mean and standard
# deviation values
		Score[,j] <- (Score[,j]-meanSx)/Sstdev

# output information on these parameters if desired by the user
		if (verbose) {
			fn1 <- sampleNames(afbatch)[idx1]
			fn2 <- sampleNames(afbatch)[idx2]
			chip <- whatcdf(fn1)
			num.probesets <- nrow(Score)
 			writeLines("S-score parameters:") 
 			writeLines(sprintf("Probearray type:      %s", chip)) 
 			writeLines(sprintf("sample1:      %s", fn1)) 
 			writeLines(sprintf("sample2:      %s", fn2))
 			writeLines(sprintf("Alpha--error coupling factor within a probeset:     %8.3f",m.alpha))
 			writeLines(sprintf("Gamma--weight of multiplicative error:     %8.3f",m.gamma))
 			writeLines(sprintf("Number of probesets:      %d",num.probesets))
	 		writeLines(" ")
 		 	writeLines(sprintf("Scaling Factor:      %8.3f     %8.3f",SF[idx1],SF[idx2]))
 			writeLines(sprintf("SDT background noise:      %8.3f     %8.3f",SDT[idx1],SDT[idx2]))
 			writeLines(sprintf("Max Intensity:      %8.3f     %8.3f",max1,max2))
 			writeLines(" ")
		}
	}

# round the S-Scores and CorrDiff to the number of digits specified by
# the user.  For the desktop version, this was 3
	if (!is.null(digits)) {
		Score <- round(Score,digits)
		CorrDiff <- round(CorrDiff,digits)
	}

	fnames <- sampleNames(afbatch)
	dimnames(Score) <- list(geneNames(afbatch),paste(fnames[compare[,1]],"vs",fnames[compare[,2]]))
	dimnames(CorrDiff) <- list(geneNames(afbatch),paste(fnames[compare[,1]],"vs",fnames[compare[,2]]))

# put the values into an ExprSet to return.  The phenoData, annotation,
# and description are the same as the AffyBatch object
	eset <- new("exprSet",
		exprs=Score,
		se.exprs=CorrDiff,
		phenoData=phenoData(afbatch),
		annotation=annotation(afbatch),
		description=description(afbatch),
		notes=notes(afbatch))
	return(eset)
}

