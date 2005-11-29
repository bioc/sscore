OneStepBiweightAlgorithm <- function(x,c,epsilon) {
#######################################################################
#
# This function computes one-step Tukey's biweight on a vector.  Note
# that this implementation follows the Affymetrix code, which is
# different from the Tukey's biweight computed by the affy package.
#
# Arguments:
#       x - vector of data on which to compute biweight value
#       c - tuning constant described in Affymetrix Statistical
#           Algorithms Description Document
#       epsilon - fuzz value to avoid division by zero described in
#                 Affymetrix Satistical Algorithms Description Document
#
# Value:
#	a numeric value of the Tukey's biweight value of the cor-
#       responding data vector
#
#######################################################################

	medianValue <- median(x)
	mad.value <- median(abs(x - medianValue))*c + epsilon
	n <- length(x)
	diff <- x - medianValue
	u <- diff / mad.value
	uSquare <- u * u
	oneMinusUSquare <- 1 - uSquare
	weightedSumNumer <- sum(ifelse(abs(u) < 1, diff * oneMinusUSquare * oneMinusUSquare, 0))
	weightedSumDenom <- sum(ifelse(abs(u) < 1, oneMinusUSquare * oneMinusUSquare, 0))
	if (weightedSumDenom != 0)
		value <- medianValue + weightedSumNumer / weightedSumDenom else
		value <- 0
	return(value)
}

#######################################################################
# The corresponding Affymetrix C++ code follows below
#######################################################################
# float OneStepBiweightAlgorithm(const vector<float> & x, float c, float epsilon)
# {
#	float medianValue = median(x);	
#	float MAD = medianAbsoluteDeviation(x) * c + epsilon;
#	int n = x.size();
#	float value=0.0f;
#	float weightedSumNumer=0.0f;
#	float weightedSumDenom=0.0f;
#
#	for (int i=0; i<n; i++)
#	{
#		float diff = x[i] - medianValue;
#		float u = diff / MAD;
#		float uSquare = u * u;
#		float oneMinusUSqaure = 1.0 - uSquare;
#		if (fabs(u) < 1.0f)
#		{
#			weightedSumNumer += diff * oneMinusUSqaure * oneMinusUSqaure;
#			weightedSumDenom += oneMinusUSqaure * oneMinusUSqaure;
#		}
#	}
#	
#	if (weightedSumDenom != 0.0f)
#	{
#		value = medianValue + weightedSumNumer / weightedSumDenom;
#	}
#
#	return value;
# }
#######################################################################

