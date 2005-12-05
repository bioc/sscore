computeOutlier <- function(afbatch,rm.mask = TRUE, rm.outliers = TRUE, rm.extra = TRUE,
  celfile.path = NULL) {
#######################################################################
#
# Arguments:
#	afbatch - an AffyBatch object containing the names of the chips
#                  to calculate outliers
#	rm.outliers,
#       rm.mask,
#       rm.extra - logical value to remove outliers and/or masked
#                  values.  If rm.extra = TRUE, it overrides the values
#                  of rm.outliers and rm.mask.  Works identically to
#                  these parameters in ReadAffy(), which it calls
#       celfile.path - character string giving the directory path to
#                      the *.CEL files being read, or NULL if the *.CEL
#                      files are in the current directory
#
# Value:
#	a matrix containing the list of outliers / masked values for
#       the given AffyBatch object.  The number of rows in the matrix
#       is equal to the number of probes for a .CEL file, and the
#       number of columns is equal to the number of chips (columns of
#       AffyBatch).  The value of each location in the matrix will be
#       TRUE if the corresponding probe is an outlier / masked value
#       and FALSE if it is not.  The probes will be arranged in the
#       same order as the intensity values, so that the outliers
#       belonging to a specific probeset can be accessed using the
#       pmindex / mmindex functions.  Note that this function assumes
#       the .CEL files are still available in the current directory.
#
# Data dictionary:
#	filenames - a list of the filenames corresponding to the
#                   AffyBatch object, in the same order as the
#                   AffyBatch columns, for reading in the .CEL file
#                   data
#	cel - an AffyBatch object containing the .CEL file data with
#             the outliers flagged
#######################################################################

# get the filenames corresponding to the AffyBatch object
	filenames <- sampleNames(afbatch)
	writeLines("Computing outliers")

# read in the .CEL files again, omitting the outliers.  This has to be
# done because the ReadAffy() function returns NAs for the outliers, so
# we cannot access the values of outliers in later calculations
# otherwise
	cel <- ReadAffy(filenames=filenames,rm.mask=rm.mask,rm.outliers=rm.outliers,rm.extra=rm.extra,celfile.path=celfile.path)

# set up the outlier matrix with TRUE/FALSE values
	outlier <- is.na(intensity(cel))
	return(outlier)
}
