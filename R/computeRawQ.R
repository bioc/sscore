computeRawQ <- function(fname,intensity,probe.index,probe.zoneID,
  bgCells,NumberZones,celfile.path = NULL) {
#######################################################################
#
# This function computes the RawQ value of a single Affymetrix GeneChip
#
# Arguments:
#       fname - character string with the filename of the GeneChip
#       intensity - vector of intensities for the GeneChip
#       probe.index - vector of indices for each probe
#       probe.zoneID - vector of zone ID numbers for each probe
#       bgCells - number of background cells for the GeneChip
#       NumberZones - number of zones on the GeneChip
#       celfile.path - character string giving the directory path to
#                      the *.CEL files being read, or NULL if the *.CEL
#                      files are in the current directory
#
# Value:
#	a numeric value corresponding to the RawQ value for the given
#       array
#
#######################################################################

# read in the stdv and npixels column of the *.CEL file, which are not
# stored in the AffyBatch object
	if (is.null(celfile.path))
		filename <- fname else
		filename <- file.path(celfile.path,fname)
	ff <- file(filename)
	open(ff,"r")
	dataline <- scan(ff,what=character(1),nmax=1,flush=TRUE,quiet=TRUE)
	while (dataline != "[INTENSITY]")
		dataline <- scan(ff,what=character(1),nmax=1,flush=TRUE,quiet=TRUE)
	dataline <- scan(ff,what=character(1),nmax=1,flush=TRUE,quiet=TRUE)
	numCells <- as.integer(unlist(strsplit(dataline,"="))[[2]])
	dataline <- scan(ff,what=character(1),nmax=1,flush=TRUE,sep="\n",quiet=TRUE)
	data <- scan(ff,what=list(NULL,NULL,NULL,stdv=0,npixels=0),nmax=numCells,flush=TRUE,quiet=TRUE)

# find the noise for each zone
	zonenoise <- tapply(probe.index,probe.zoneID,computeZonenoise,intensity,data$stdv,data$npixels,bgCells)

# rawQ is the average noise over all zones
	rawQ <- sum(zonenoise) / NumberZones
	return(rawQ)
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

