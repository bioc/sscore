computeZonenoise <- function(index,intensity,stdv,npixels,bgCells) {
#######################################################################
#
# This function computes the noise (average standard error) of the
# probe intensities for a single Affymetrix GeneChip
#
# Arguments:
#       index - vector of indices for probes in the given zone
#       intensity - vector of intensities for the GeneChip
#       stdv - vector of standard deviations for the GeneChip
#       npixels - vector containing number of pixels for each probe of
#                 the GeneChip
#       bgCells - number of background cells on the GeneChip
#
# Value:
#	an object of class ExprSet, with the exprs slot containing the
#       S-score values and the se.exprs slot containing the CorrDiff
#       values
#
#######################################################################

	cellIndex <- index[order(intensity[index])[1:bgCells]]
	bgSum <- sum(stdv[cellIndex] / sqrt(npixels[cellIndex]))
	return(bgSum / bgCells)
}

#######################################################################
# The corresponding Affymetrix C++ code follows below
#######################################################################
# float CRawQ::ComputeRawQ(CCELFileData &cell, CCDFFileData &cdf)
# {
#	typedef struct {
#		float intensity;
#		float stdev;
#		int   pixel;
#	} CellStatisticsType;
#
#
#	// Store the number of background cells.
#	CCDFFileHeader &cdfHeader = cdf.GetHeader();
#	float numer = cdfHeader.GetCols() * cdfHeader.GetRows() * m_PercentBGCells;
#	float denom = m_VertZones * m_HorzZones * 100.0f;
#	int bgCells = (int) (numer / denom);
#
#
#	// Determine the number of remaining cells in the vertical direction.
#	int CellsRemaining = cdfHeader.GetCols() % m_VertZones;
#	int zonex;
#	int zoney;
#	if(CellsRemaining != 0)
#		zonex = (cdfHeader.GetCols() + 
#				m_VertZones - CellsRemaining) /
#				m_VertZones;
#	else
#		zonex = cdfHeader.GetCols() / m_VertZones;
#
#	// Determine the number of remaining cells in the horizontal direction.
#	CellsRemaining = cdfHeader.GetRows() % m_HorzZones;
#	if(CellsRemaining != 0)
#		zoney = (cdfHeader.GetRows() +
#				m_HorzZones-CellsRemaining) /
#				m_HorzZones;
#	else
#		zoney = cdfHeader.GetRows() / m_HorzZones;
#
#	// Ensure that there are a match and mismatch cell in the same zone.
#	zoney += zoney % EXPRESSION_ATOMS_PER_CELL;
#
#
#	// Determine the total number of zones.
#	int NumberZones = m_VertZones * m_HorzZones;
#
#	// Allocate memory for storing background data.
#	float *zonebg = new float[NumberZones];
#	float *zonenoise = new float[NumberZones];
#	CellStatisticsType *bgN = new CellStatisticsType[NumberZones*bgCells];
#
#	// Get number of units.
#	int NumUnits = cdfHeader.GetNumProbeSets();
#
#	// Determine the total number of atoms in the chip.
#	CCDFProbeSetInformation unit;
#	int totalCells=0;
#	for (int iUnit=0; iUnit<NumUnits; ++iUnit)
#	{
#		if (cdf.GetProbeSetType(iUnit) == ExpressionProbeSetType)
#		{
#			cdf.GetProbeSetInformation(iUnit, unit);
#			totalCells += unit.GetNumCells();
#		}
#	}
#
#	// Allocate space for all atoms intensities and ID's.
#	int *entryIndex = new int[totalCells];
#	int *intenszid = new int[totalCells];
#
#	// Clear arrays.
#	memset(intenszid, 0, sizeof(int)*totalCells);
#	memset(entryIndex, 0, sizeof(int)*totalCells);
#
#	// Loop over all units to determine the zone ID's and intensities.
#	int iInten=0;
#	for (int iUnit=0; iUnit<NumUnits; ++iUnit)
#	{
#		// Only process expression units.
#		if (cdf.GetProbeSetType(iUnit) != ExpressionProbeSetType)
#			continue;
#
#		// Get the PM and MM intensity objects and their zone ids.
#		cdf.GetProbeSetInformation(iUnit, unit);
#		CCDFProbeGroupInformation group;
#		int nGroups = unit.GetNumGroups();
#		for (int iGroup=0; iGroup<nGroups; iGroup++)
#		{
#			unit.GetGroupInformation(iGroup, group);
#			int nCells = group.GetNumCells();
#			CCDFProbeInformation probe;
#			for (int iCell=0; iCell<nCells; iCell++)
#			{
#				group.GetCell(iCell, probe);
#				entryIndex[iInten] = cell.XYToIndex(probe.GetX(), probe.GetY());
#				intenszid[iInten] = DetermineZone(probe.GetX(), probe.GetY(), zonex, zoney, m_VertZones);
#				++iInten;
#			}
#		}
#	}
#
#	// compute background for each zone
#	for (int iZone=0; iZone<NumberZones; iZone++)
#	{
#		// Initialize the background.
#		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
#			bgN[bgcnt+(iZone*bgCells)].intensity = LARGE_FLOAT_NUMBER;
#
#		// find the lowest N intensities in each zone
#		for (int iInten = 0; iInten < totalCells; iInten++)
#		{
#			// Only process those intensities in the current zone.
#			if(intenszid[iInten] == iZone)
#			{
#				int index_cnt;
#				int index_max;
#				for (index_cnt=1, index_max=0;
#					 index_cnt < bgCells; index_cnt++)
#				{
#					if(bgN[index_cnt+(iZone*bgCells)].intensity > bgN[index_max+(iZone*bgCells)].intensity)
#						index_max = index_cnt;
#				}
#
#
#				// Store the low intensity.
#				float intensity = min(bgN[index_max+(iZone*bgCells)].intensity, cell.GetIntensity(entryIndex[iInten]));
#				if (intensity != bgN[index_max+(iZone*bgCells)].intensity)
#				{
#					bgN[index_max+(iZone*bgCells)].intensity = cell.GetIntensity(entryIndex[iInten]);
#					bgN[index_max+(iZone*bgCells)].pixel = cell.GetPixels(entryIndex[iInten]);
#					bgN[index_max+(iZone*bgCells)].stdev = cell.GetStdv(entryIndex[iInten]);
#				}
#			}
#		}
#
#		// compute the average
#		float bgSum = 0.0f;
#		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
#			bgSum += bgN[bgcnt+(iZone*bgCells)].intensity;
#		zonebg[iZone] = bgSum / bgCells;
#
#		// Compute the noise.
#		bgSum = 0.0f;
#		for(int bgcnt = 0; bgcnt < bgCells; bgcnt++)
#			bgSum += bgN[bgcnt+(iZone*bgCells)].stdev / (float)(sqrt((float)bgN[bgcnt+(iZone*bgCells)].pixel));
#		zonenoise[iZone] = bgSum / bgCells;
#	}
#
#	// Compute the average noise.
#	float avgNoise=0.0f;
#	for (int iZone=0; iZone<NumberZones; iZone++)
#		avgNoise+=zonenoise[iZone];
#	float rawQ = avgNoise/NumberZones;
#
#	//Clean up
#	delete [] entryIndex;
#	delete [] intenszid;
#	delete [] bgN;
#	delete [] zonebg;
#	delete [] zonenoise;
#
#	// Return the RawQ value
# 	return rawQ;
# }
#######################################################################

