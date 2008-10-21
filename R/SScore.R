SScore <- function(afbatch = stop("No CEL files specified"), classlabel =
 c(0,1), SF = NULL, SDT = NULL, rm.outliers = TRUE, rm.mask = TRUE,
 rm.extra = NULL, digits = NULL, verbose = FALSE,celfile.path = NULL, 
 celfile.names=NULL) {


#######################################################################
#
# This function computes the S-Score values for each probeset in an
# AffyBatch object.  S-Scores are estimated for one pair of GeneChips.
#
# Arguments:
#	afbatch - an AffyBatch object containing the data on the chips
#                 to be compared
#	classlabel - a numeric vector specifying the classes to be compared.
#                Its length is equal to the number of chips in afbatch,
#                with each element being a 0 or 1 to denote class
#                membership of that chip.  S-Scores will be generated
#                comparing all class 0 chips against all class 1 chips.
#                The default value of c(0,1) will give a 2-chip 
#                comparisons compatible with previous releases of the
#                function.
#	SF     -  a vector of SF values (the Scale Factor, which is
#                 part of the general output from the Affymetrix GCOS
#                 software), one for each chip.  If SF = NULL, these
#                 will be calculated internally
#   SDT    -  a vector of SDT values (the Standard Difference
#                 Threshold, which can be derived from the general 
#                 output from the Affymetrix GCOS software using the
#                 formula SDT = 4 * RawQ * Scale Factor), one for each
#                 chip.  If SDT = NULL, these will be calculated
#                 internally
#	rm.outliers,
#   rm.mask,
#   rm.extra - logical value to remove outliers and/or masked 
#                  values.  If rm.extra = TRUE, it overrides the values
#                  of rm.outliers and rm.mask.  Works identically to
#                  these parameters in ReadAffy(), which it calls
#	digits - the number of significant digits to retain when
#                returning S-Score values
#	verbose  - logical value indicating whether additional informa-
#                  tion on the calculations should be displayed
#   celfile.path - character string giving the directory path to
#                      the *.CEL files being read, or NULL if the *.CEL
#                      files are in the current directory
#
# Value:
#	an object of class ExpressionSet, with the exprs slot containing the
#       S-score values and the se.exprs slot containing the CorrDiff
#       values
#
# Data dictionary:
#	sfsdtlist - internal list of calculated SF and SDT values, if
#               needed	
#	outlier - matrix of logicals identifying the outlier / masked
#             values.  It is arranged in the same format as the
#             intensities (# of rows equal to number of probes, #
#             of columns equal to number of chips) with each entry
#             being TRUE if the probe for that chip is an outlier/
#		      masked value or FALSE if it is not.
#	m.gamma - value of gamma in the S-score equation (weight of
#             multiplicative error)
#	m.alpha - value of alpha in the S-Score equation (coupling
#             factor among probe pairs in probeset)
#	probenames - list of the probe set names for the set of chips
#	pmidx,
#	mmidx - list of the PM and MM indices, respectively, for the
#           set of chips
#	Score - matrix containing the S-Score values
#	CorrDiff - matrix containing the CorrDiff values
#	idx1,
#	idx2 - indices into afbatch for the two chips currently being
#          compared
#	intens1,
#	intens2 - intensity values for the two chips being compared
#	max1,
#	max2 - maximum intensity (for PM and MM combined) for the two
#          chips being compared
#	PM1,
#	PM2 - PM intensity values of the current probeset for the two
#         chips being compared
#	MM1,
#	MM2 - MM intensity values of the current probeset for the two
#         chips being compared
#	index - indices for the intensity values of the current probe-
#           set that are being used in the S-Score calculation
#           (i.e., not saturated and if applicable, not outlier /
#		    masked values)
#	outlier1,
#	outlier2 - logical vector indicating outlier / masked values of
#              the current probeset for the two chips being
#              compared
#	N - number of probes of the current probeset that are not
#       saturated (as well as not outlier / masked values).
#	diff1,
#	diff2 - PM-MM differences for the two chips being compared
#	f.err - "raw" S-Score value (without N or alpha)
#	Sx - sum of x  (chip1 PM-MM differences or chip pair S-Scores)
#        values
#	Sxx - sum of x-squared (variance for chip 1 Pm-MM differences
#         or chip pair S-Scores) values
#	Sy - sum of y (chip 2 PM-MM differences) values
#	Syy - sum of y-squared (variance for chip 2 PM-MM differences)
#         values
#	Sxy - sum of covariance values between chips 1 and 2 PM-MM
#         differences
#	x - temporary variable for storing S-Score values in
#       calculations
#	num - number of probe pairs in trimmed probeset
#	meanSx - mean of x (chip pair S-Scores) values
#	Sstdev - standard deviation of x (chip pair S-Scores) values
#	fn1,
#	fn2 - filenames of the two chips being compared
#	chip - name of CDF / type of chips being compared
#	i - counter for probeset number in loops
#
########################################################################
# Version history
#
# Note: old code is preceded by "###" to identify removals
########################################################################
# check the number of CEL files
###	if (length(afbatch) < 2)
###		stop("Must have two chips for comparisons")
###	if (length(afbatch) > 2)
###		warning("More than two chips listed for comparison.  Only the first two will be used")
	bad.label <- is.na(match(classlabel,c(0,1)))
	if (any(bad.label)) 
		stop("Two classes, labeled 0 and 1, are required for comparisons")

	if (is.null(celfile.names))
		fname <- sampleNames(afbatch) else
		fname <- celfile.names
	if (length(fname) != length(afbatch))
		stop("Number of filenames does not match number of samples")
		
	outlier <- matrix(data=FALSE,nrow=nrow(intensity(afbatch)),ncol=ncol(intensity(afbatch)))
	if (!is.null(rm.extra)) 
		rm.outliers <- rm.mask <- rm.extra else
		rm.extra <- FALSE
	if (is.null(SF) | is.null(SDT) | rm.outliers | rm.mask | rm.extra) {
		stdvs <- pixels <- NULL
		for (i in 1:length(fname)) {
###			celdata <- readCel(fname[i],readHeader=FALSE,readIntensities=FALSE,readStdvs=TRUE,readPixels=TRUE)
			if (is.null(celfile.path))
				filename <- fname[i] else
				filename <- file.path(celfile.path,fname[i])
			celdata <- read.celfile(filename)
			if (is.null(SF) | is.null(SDT)) 
###				stdvs <- cbind(stdvs,celdata$stdvs)
				stdvs <- cbind(stdvs,celdata$INTENSITY$STDEV)
###				pixels <- cbind(pixels,celdata$pixels)
				pixels <- cbind(pixels,celdata$INTENSITY$NPIXELS)
			if (rm.outliers | rm.extra) {
				writeLines("Computing outliers")
				outlier.xy <- celdata$OUTLIERS
				outlier[xy2indices(outlier.xy[,1],outlier.xy[,2],abatch=afbatch),i] <- TRUE
			}
			if (rm.mask | rm.extra) {
				outlier.xy <- celdata$OUTLIERS
				outlier[xy2indices(outlier.xy[,1],outlier.xy[,2],abatch=afbatch),i] <- TRUE
			}
		}
		if (is.null(SF) | is.null(SDT)) {
			writeLines("Computing SF and SDT")
			sfsdtlist <- computeAffxSFandSDT(afbatch,stdvs,pixels,digits=3)
		}
		if (is.null(SF))
			SF <- sfsdtlist$SF
		if (is.null(SDT))
			SDT <- sfsdtlist$SDT
	}

# calculate SF and SDT, if not specified by the user, as these must
# always be available
	if (any(SF <= 0)) 
		stop("SF values must be positive")
###	if (length(SF) != 2)
###		stop("Must be two SF values, one for each CEL file")
	if (length(SF) != length(afbatch))
		stop("Number of SF values does not match number of samples")

	if (any(SDT <= 0)) 
		stop("SDT values must be positive")
###	if (length(SDT) != 2)
###		stop("Must be two SDT values, one for each CEL file")
	if (length(SDT) != length(afbatch))
		stop("Number of SDT values does not match number of samples")

# calculate the outlier matrix, if excluding outliers / masked values
#	if (any(rm.outliers,rm.mask,rm.extra))
#		outlier <- computeOutlier(afbatch,rm.outliers=rm.outliers,rm.mask=rm.mask,rm.extra=rm.extra,celfile.path=celfile.path)

# initialize variables needed in later calculations.  For the PM/MM
# index, intens1/intens2, and max1/max2 values, these are calculated
# outside of the loop, since they are constant for a given CDF or chip.
	m.gamma <- 0.1
	probenames <- featureNames(afbatch)
	pmidx <- pmindex(afbatch)
	mmidx <- mmindex(afbatch)
###	intens1 <- intensity(afbatch[,1])*SF[1]
###	intens2 <- intensity(afbatch[,2])*SF[2]
	intens1 <- t(t(intensity(afbatch[,classlabel==0]))*SF[classlabel==0])
	intens2 <- t(t(intensity(afbatch[,classlabel==1]))*SF[classlabel==1])
	Score <- CorrDiff <- rep(0.0,length(probenames))
	writeLines("Computing S-score values")
	max1 <- apply(rbind(pm(afbatch[,classlabel==0]),mm(afbatch[,classlabel==0])),2,max)*SF[classlabel==0]
	max2 <- apply(rbind(pm(afbatch[,classlabel==1]),mm(afbatch[,classlabel==1])),2,max)*SF[classlabel==1]

# this loops through each of the probesets on a given pair of chips
###	for (i in 1:length(probenames)) {
	for (i in 1:length(pmidx)) {

# get the PM and MM values, as well as minimum intensity values, for
# the given probeset on each chip of the pair
###		PM1 <- intens1[pmidx[[i]]]
###		MM1 <- intens1[mmidx[[i]]]
###		PM2 <- intens2[pmidx[[i]]]
###		MM2 <- intens2[mmidx[[i]]]
		PM1 <- intens1[pmidx[[i]],,drop=FALSE]
		MM1 <- intens1[mmidx[[i]],,drop=FALSE]
		PM2 <- intens2[pmidx[[i]],,drop=FALSE]
		MM2 <- intens2[mmidx[[i]],,drop=FALSE]
###		min1 <- min(PM1,MM1)
###		min2 <- min(PM2,MM2)
		min1 <- apply(rbind(PM1,MM1),2,min)
		min2 <- apply(rbind(PM2,MM2),2,min)
 	
# adjust each of the PM and MM intensities relative to the minimum
# values
###		PM1 <- PM1 - min1 
###		PM2 <- PM2 - min2 
###		MM1 <- MM1 - min1 
###		MM2 <- MM2 - min2
		PM1 <- t(t(PM1) - min1) 
		PM2 <- t(t(PM2) - min2) 
		MM1 <- t(t(MM1) - min1) 
		MM2 <- t(t(MM2) - min2)

# find the index of the probe pairs of the probeset to use in
# calculations.  A probe pair is used if it is not "saturated" (i.e.,
# the intensity is less than the maximum - minimum) and if it is not
# identified as an outlier / masked value in the .CEL file
		index <- cbind((PM1<max1-min1),(PM2<max2-min2),(MM1<max1-min1),(MM2<max2-min2))
		index <- apply(index,1,any)
		if (any(rm.outliers,rm.mask,rm.extra)) {
###			outlier1 <- outlier[pmidx[[i]],1] | outlier[mmidx[[i]],1]
###			outlier2 <- outlier[pmidx[[i]],2] | outlier[mmidx[[i]],2] 
###			index <- index & (!outlier1) & (!outlier2)
			outlier1 <- outlier[pmidx[[i]],classlabel==0,drop=FALSE] | outlier[mmidx[[i]],classlabel==0,drop=FALSE]
			outlier2 <- outlier[pmidx[[i]],classlabel==1,drop=FALSE] | outlier[mmidx[[i]],classlabel==1,drop=FALSE] 
			index <- index & (!apply(outlier1,1,any)) & (!apply(outlier2,1,any))
		}
		N <- sum(as.integer(index))

# find the PM-MM differences for the probeset on each of the two chips in the pair
###		diff1 <- (PM1-MM1)[index] 
###		diff2 <- (PM2-MM2)[index] 
		diff1 <- (PM1-MM1)[index,,drop=FALSE] 
		diff2 <- (PM2-MM2)[index,,drop=FALSE] 

# this is the "raw" S-Score value (without N or alpha in calculation) -
# compare to S-Score calculation in J Mol Biol paper
		N1 <- sum(classlabel==0)
		N2 <- sum(classlabel==1)
		f.err <- (apply(diff1,1,sum)/N1-apply(diff2,1,sum)/N2)/
		  sqrt(m.gamma*m.gamma*(apply(diff1*diff1,1,sum)/(N1*N1)+apply(diff2*diff2,1,sum)/(N2*N2))+
		  sum(SDT[classlabel==0]*SDT[classlabel==0])/(N1*N1)+sum(SDT[classlabel==1]*SDT[classlabel==1])/(N2*N2)) 

# threshold or impute outlying f.err values.  The cutoff of 15 was
# arbitrarily decided in the original version; how it was determined
# is unknown
		f.err[f.err > 15.0] <- 15.0 
		f.err[f.err < -15.0] <- -15.0 
 
		Score[i] <- sum(f.err) 
 
# estimate the variance / covariance values, for calculating the
# CorrDiff
###		Sxx <- sum(diff1*diff1) 
###		Syy <- sum(diff2*diff2) 
###		Sxy <- sum(diff1*diff2) 
		Sxx <- sum(mean(diff1)^2) 
		Syy <- sum(mean(diff2)^2) 
		Sxy <- sum(mean(diff1)*mean(diff2)) 
 
		Sx <- 0 
		Sy <- 0 

# transform the S-Score estimate by dividing by a function of the
# number of probes in the probeset
		if (N > 0) 
			Score[i] <- Score[i] / sqrt(N) else Score[i] <- 0 
 
# calculate the CorrDiff.  CorrDiffs below the threshold of 1e-3 are
# imputed to be 0, which was also arbitrarily decided in the first
# version
		if (N>2 && ((Sxx-Sx*Sx/N)*(Syy-Sy*Sy/N) > 1.e-3)) 
			CorrDiff[i] <- (Sxy-Sx*Sy/N)/sqrt((Sxx-Sx*Sx/N)*(Syy-Sy*Sy/N)) else CorrDiff[i] <- 0.0 
	}

# now renormalize the S-Score values, which gives alpha
	writeLines("Renormalizing S-scores")
	x <- Score
	Sx <- sum(x)
	Sxx <- sum(x*x)

# calculate the mean and standard deviation of the entire set of
# S-Scores
	Sstdev <- sqrt((Sxx-Sx*Sx/length(Score))/length(Score))
	meanSx <- Sx/length(Score)

# find the trimmed S-score values, using a cutoff of those S-Scores
# within 3 standard deviations of the mean
	x <- Score-meanSx;
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
	Score <- (Score-meanSx)/Sstdev

###	fn1 <- fname[1]
###	fn2 <- fname[2]
	fn1 <- paste(fname[classlabel==0],collapse=" ")
	fn2 <- paste(fname[classlabel==1],collapse=" ")

# output information on these parameters if desired by the user
	if (verbose) {
		Chip <- cdfName(afbatch)
		num.probesets <- length(Score)
 		writeLines("S-score data. Parameter section:") 
 		writeLines(sprintf("Probearray type:      %s", Chip)) 
 		writeLines(sprintf("sample1:      %s", fn1)) 
 		writeLines(sprintf("sample2:      %s", fn2))
 		writeLines(sprintf("Alpha--error coupling factor within a probeset:     %8.3f",m.alpha))
 		writeLines(sprintf("Gamma--weight of multiplicative error:     %8.3f",m.gamma))
 		writeLines(sprintf("Number of Probesets:     %i",num.probesets))
		writeLines(" ")
 		writeLines("Scaling Factor:")
		printSF <- formatC(SF[classlabel==0],digits=3,width=8,format="f")
		writeLines(sprintf("  sample1 (class label 0):      %s",paste(printSF,collapse=" ")))
		printSF <- formatC(SF[classlabel==1],digits=3,width=8,format="f")
		writeLines(sprintf("  sample2 (class label 1):      %s",paste(printSF,collapse=" ")))
 		writeLines("SDT background noise:")
 		printSDT <- formatC(SDT[classlabel==0],digits=3,width=8,format="f")
 		writeLines(sprintf("  sample1 (class label 0):      %s",paste(printSDT,collapse=" ")))
 		printSDT <- formatC(SDT[classlabel==1],digits=3,width=8,format="f")
 		writeLines(sprintf("  sample2 (class label 1):      %s",paste(printSDT,collapse=" ")))
 		writeLines("Max Intensity:")
 		printMax <- formatC(max1,digits=3,width=8,format="f")
 		writeLines(sprintf("  sample1 (class label 0):      %s",paste(printMax,collapse=" ")))
 		printMax <- formatC(max2,digits=3,width=8,format="f")
 		writeLines(sprintf("  sample2 (class label 1):      %s",paste(printMax,collapse=" ")))
 		writeLines(" ")
	}

# round the S-Scores and CorrDiff to the number of digits specified by
# the user.  For the desktop version, this was 3
	if (!is.null(digits)) {
		Score <- round(Score,digits)
		CorrDiff <- round(CorrDiff,digits)
	}

	Score <- as.matrix(Score)
	CorrDiff <- as.matrix(CorrDiff)
	rownames(Score) <- rownames(CorrDiff) <- featureNames(afbatch)
#	colnames(Score) <- colnames(CorrDiff) <- paste(fn1,"vs",fn2)
	colnames(Score) <- colnames(CorrDiff) <- "Class 0 vs 1"
	comparison <- 1
	Score.pData <- data.frame(comparison,row.names="Class 0 vs 1")
#	ScorePheno <- new('phenoData', pData=Score.pData,
#	  varLabels=list(sample = "arbitrary numbering"))
	Score.Metadata <- data.frame(labelDescription = "arbitrary numbering", 
	  row.names = "comparison")
	ScorePheno <- new("AnnotatedDataFrame", data=Score.pData, varMetadata =
	  Score.Metadata)
# put the values into an ExprSet to return.  The phenoData, annotation,
# and description are the same as the AffyBatch object
	eset <- new("ExpressionSet",
		exprs=Score,
#		se.exprs=CorrDiff,
		phenoData=ScorePheno,
		annotation=annotation(afbatch))
#		description=description(afbatch),
#		notes=notes(afbatch))
	return(eset)
}
