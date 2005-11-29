trimMean <- function(vec,p1,p2) {
#######################################################################
#
# This function computes the trimmed mean for a vector.  Note that this
# implementation follows the Affymetrix code, which gives different
# results than the standard R function mean().
#
# Arguments:
#       vec - vector of values for computing trimmed mean
#       p1,
#       p2 - lower and upper percentage for trimming, expressed as a
#            decimal fraction (not whole number)
#
# Value:
#	a numeric value representing the trimmed mean for the given
#       data vector
#
#######################################################################

	whole <- vec
	total <- length(vec)
	if (total==0) return(0)

	whole <- sort(whole)

	dG1 <- total * p1 + 1
	dG2 <- total * (1 - p2) + 1
	g1 <- floor(dG1)
	g2 <- floor(dG2)
	r1 <- dG1 - g1
	r2 <- dG2 - g2
	last <- total - g2 + 1
	if (last <= 0) last <- 0

	sum <- (1 - r1) * whole[g1] + (1 - r2) * whole[last]
	sum <- sum + sum(whole[(g1+1):(last-1)])

	subtotal <- last - g1 - 1
	if (subtotal <= 0) subtotal <- 0
	subtotal <- subtotal + 2 - r1 - r2
	return(sum / subtotal)
}

#######################################################################
# The corresponding Affymetrix C++ code follows below
#######################################################################
# template <class T> double trimMean(const vector<T> & vec, const double p1, const double p2) 
# {
#	vector<T> whole = vec;
#	int total = whole.size();
#	if (total == 0)
#		return 0.0f;
#
#	sort(whole.begin(), whole.end());
#
#	double dG1 = total * p1;
#	double dG2 = total * (1.0 - p2);
#	int g1 = floor(dG1);
#	int g2 = floor(dG2);
#	double r1 = dG1 - g1;
#	double r2 = dG2 - g2;
#	int last = total - g2 - 1;
#	if (last <= 0.0f) { // it is theoretically impossible for last < 0, but
#		last = 0; // we add the code here to guarantee proper bounds even if there is any numerical unstability
#	}
#	double sum = (1.0f - r1) * whole[g1] + (1.0f - r2) * whole[last];
#	for (int i = g1 + 1; i < last; ++i) {
#		sum += whole[i];
#	}
#	double subtotal = last - g1 -1;
#	if (subtotal <= 0.0f) {
#		subtotal = 0.0;
#	}
#	subtotal += 2.0 - r1 - r2;
#	return sum / subtotal;
# }
#######################################################################
