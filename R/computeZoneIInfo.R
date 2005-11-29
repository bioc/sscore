computeZoneIInfo <- function(ZoneInfo,NumberBGCells) {
#######################################################################
#
# This function computes the background and noise for a given zone of a
# single Affymetrix GeneChip
#
# Arguments:
#       ZoneInfo - vector of intensities in a given zone
#       NumberBGCells - number of background cells for the GeneChip
#
# Value:
#       a list containing the components background (background value
#       for the given zone)  and noise (noise value for the given zone)
#
#######################################################################

	lowBG <- 0
	highBG <- NumberBGCells / 100
	if (length(ZoneInfo) == 0) 
		ZonesInfo.background <- ZonesInfo.noise <- 0 else {
			sort.list <- sort(ZoneInfo)
			n1 <- 1
			n2 <- floor(length(ZoneInfo) * highBG)
			subtotal <- n2
			if (subtotal > 2) {
				sum <- sum(sort.list[n1:n2])
				tMean <- sum / subtotal
				sum <- sum((sort.list[n1:n2] - tMean)^2)
				ZonesInfo.background <- tMean
				ZonesInfo.noise <- sqrt(sum / (subtotal-1))
			} else {
				ZonesInfo.background <- sort.list[n1]
				ZonesInfo.noise <- 0
			}	
		}
	ZoneIInfo <- c(ZonesInfo.background,ZonesInfo.noise)
	names(ZoneIInfo) <- c("background","noise")
	return(ZoneIInfo)
}

#######################################################################
# The corresponding Affymetrix C++ code follows below
#######################################################################
# template <class T> FloatPair trimMeanAndStd(vector<T> & v, const double p1, const double p2) {
#	int total = v.size();
#	FloatPair fp(0, 0);
#	if (total > 0) {
#		sort(v.begin(), v.end());
#		int n1 = 0;
#		int n2 = floor(total * p2);
#		double subtotal = n2;
#		if (subtotal > 2) {
#			double sum = 0;
#			int i;
#			for (i = n1; i < n2; ++i) {
#				sum += v[i];
#			}
#			double tMean = sum / subtotal;
#			sum = 0;
#			for (i = n1; i < n2; ++i) {
#				sum += pow(v[i] - tMean, 2);
#			}
#			fp.value1 = (float) tMean;
#			fp.value2 = (float) sqrt(sum / (subtotal - 1));
#		}
#		else if (subtotal == 1) {
#			fp.value1 = v[n1];
#		}
#	}
#	return fp;
# }
#######################################################################
# void CExpressionAlgorithmImplementation::ComputeScaledAdjustedIntensity(
#											  CCELFileData *pCell,
#											  vector<vector<float> > & PM,
#											  vector<vector<float> > & MM,
#											  vector<vector<bool> > &UseAtom,
#											  vector<vector<FloatPair> > & BG,
#											  vector<vector<FloatPair> > & Noise,
#											  AllZonesInfoType & ZonesInfo)
# {
#	int iUnit;
#	int CellsRemaining;
#	int zonex;
#	int zoney;
#	int NumberZones;
#	int iInten=0;
#	string probeSetName;
#	FILE *ff;
#
#	int *NumberCellsPerZone=NULL;
#
#	// Determine the number of remaining cells in the vertical direction.
#	CellsRemaining = m_Cdf.GetHeader().GetCols() % (int)m_Params.NumberVertZones;
#	if(CellsRemaining != 0)
#		zonex = (m_Cdf.GetHeader().GetCols() + 
#				(int)m_Params.NumberVertZones - CellsRemaining) /
#				(int)m_Params.NumberVertZones;
#	else
#		zonex = m_Cdf.GetHeader().GetCols() / (int)m_Params.NumberVertZones;
#
#	// Determine the number of remaining cells in the horizontal direction.
#	CellsRemaining = m_Cdf.GetHeader().GetRows() % (int)m_Params.NumberHorZones;
#	if(CellsRemaining != 0)
#		zoney = (m_Cdf.GetHeader().GetRows() +
#				(int)m_Params.NumberHorZones-CellsRemaining) /
#				(int)m_Params.NumberHorZones;
#	else
#		zoney = m_Cdf.GetHeader().GetRows() / (int)m_Params.NumberHorZones;
#
#	// Ensure that there are a match and mismatch cell in the same zone.
#	zoney += zoney % 2; //EXPRESSION_ATOMS_PER_CELL;
#
#	// Determine the total number of zones.
#	NumberZones = (int)m_Params.NumberVertZones * (int)m_Params.NumberHorZones;
#
#	// Get number of units.
#	int NumUnits = m_Cdf.GetHeader().GetNumProbeSets();
#
#	// Allocate space for all atoms intensities and ID's.
#	NumberCellsPerZone = new int[NumberZones];
#
#	// Clear arrays.
#	memset(NumberCellsPerZone, 0, sizeof(int)*NumberZones);
#
#	// Loop over all units to determine the zone ID's and intensities.
#	vector<vector<float> > ZoneCells(NumberZones);
#	CCDFProbeSetInformation unit;
#	CCDFProbeGroupInformation blk;
#	CCDFProbeInformation cell;
#	bool bMasked;
#	for (iUnit=0; iUnit<NumUnits; ++iUnit)
#	{
#		m_Cdf.GetProbeSetInformation(iUnit, unit);
#
#		// Only process expression units.
#		if (unit.GetProbeSetType() == ExpressionProbeSetType)
#		{
#			unit.GetGroupInformation(0, blk);
#			int numCells = blk.GetNumCells();
#
#			// Loop over the atoms in the unit
#			for (int iCell=0; iCell<numCells; iCell++)
#			{
#				probeSetName = m_Cdf.GetProbeSetName(iUnit);
#				blk.GetCell(iCell, cell);
#				bMasked = pCell->IsMasked(cell.GetX(), cell.GetY());
#				if (bMasked == false)
#				{
#					int nZone,Zonex,Zoney;
#					Zonex = cell.GetX();
#					Zoney = cell.GetY();
#					nZone = DetermineZone(cell.GetX(), cell.GetY(), zonex, zoney);
#					ZoneCells[nZone].resize(ZoneCells[nZone].size() + 1);
#					ZoneCells[nZone][NumberCellsPerZone[nZone]] =
#						pCell->GetIntensity(cell.GetX(), cell.GetY());
#
#					if (nZone >= 0 && nZone < NumberZones)
#						NumberCellsPerZone[nZone]++;
#				}
#			}
#		}
#	}
#
#	// Allocate zones, set smooth factor and set num zones
#	ZonesInfo.pZones = new ZoneInfo[NumberZones];
#	ZonesInfo.number_zones = NumberZones;
#	ZonesInfo.smooth_factor = m_Params.SmoothFactorBG;
#
#	// compute background for each zone
#	for (int iZone=0; iZone<NumberZones; iZone++)
#	{
#		// Compute the center coordinates of each zone.
#		// (x1,y1) is the upper left corner
#		// (x2,y2) is the lower right corner
#		float x1 = ((int) (iZone % (int)m_Params.NumberVertZones)) * zonex;
#		float y1 = ((int) (iZone / (int)m_Params.NumberVertZones)) * zoney;
#		float x2 = x1 + zonex;
#		float y2 = y1 + zoney;
#		ZonesInfo.pZones[iZone].center.x = (x1 + x2) / 2;
#		ZonesInfo.pZones[iZone].center.y = (y1 + y2) / 2;
#
#		int iCell=0;
#		int numCell = NumberCellsPerZone[iZone];
#		ZonesInfo.pZones[iZone].numCell = numCell;
#
#		vector<float> zoneI(numCell);
#		vector<int> rank(numCell);
#
#		for (int i=0; i<numCell; i++)
#		{
#			float inten = ZoneCells[iZone][i];
#			zoneI[iCell] = ModifyIntensitySlightly(inten);
#			iCell++;
#		}
#		float lowBG = 0.0f;
#		float highBG = m_Params.NumberBGCells / 100.0f;
#		FloatPair fp = trimMeanAndStd(zoneI, lowBG, highBG);
#	double zoneT;
#	zoneT = trimMean(zoneI,0.02,0.98);
#
#		ZonesInfo.pZones[iZone].background = fp.value1;
#		ZonesInfo.pZones[iZone].noise = fp.value2;
#	}
#	// End of computing background intensity and noise for each zone.
#	// Carried zones and NumberZones as the required information.
#
#	// Compute b(x,y), n(x,y), SA(x,y) which was stored in PM[i][j] and MM[i][j]
#	float smoothF = m_Params.SmoothFactorBG;
#
#	CCDFProbeInformation pmcell;
#	CCDFProbeInformation mmcell;
#	ff = fopen("Intermediate.txt","w");
#	for (iUnit=0; iUnit<NumUnits; ++iUnit)
#	{
#		m_Cdf.GetProbeSetInformation(iUnit, unit);
#
#		// Only process expression units.
#		if (unit.GetProbeSetType() == ExpressionProbeSetType)
#		{
#			unit.GetGroupInformation(0, blk);
#			int numCells = blk.GetNumCells();
#			int nAtoms = blk.GetNumLists();
#
#			PM[iUnit].resize(nAtoms, 0.0);
#			MM[iUnit].resize(nAtoms, 0.0);
#			UseAtom[iUnit].resize(nAtoms, true);
#
#			BG[iUnit].resize(nAtoms);
#			Noise[iUnit].resize(nAtoms);
#
#			// Loop over the atoms in the unit
#			int atomCount = 0;
#			for (int iCell=0, iAtom=0; iCell<numCells; iCell+=2, ++iAtom)
#			{
#				blk.GetCell(iCell, pmcell);
#				blk.GetCell(iCell+1, mmcell);
#				bMasked = (pCell->IsMasked(pmcell.GetX(), pmcell.GetY()) == true) ||
#							  (pCell->IsMasked(mmcell.GetX(), mmcell.GetY()) == true);
#				if (bMasked == true)
#				{
#						UseAtom[iUnit][iAtom] = false;
#				} 
#
#				// Set the background (with smoothing adjustment) for matchcell.
#				Coordinate Cellxy;
#				Cellxy.x = pmcell.GetX();
#				Cellxy.y = pmcell.GetY();
#				float WeightedSumBg = 0.0f;
#				float WeightedSumNoise = 0.0f;
#				float WeightedSumDenom = 0.0f;
#				float background = 0.0f;
#				float noise = 0.0f;
#				int k;
#				for (k = 0; k < NumberZones; k++)
#				{
#					WeightedSumBg    += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].background;
#					WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].noise;
#					WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF); 
#					probeSetName = m_Cdf.GetProbeSetName(iUnit);
#				}
#				probeSetName = m_Cdf.GetProbeSetName(iUnit);
#				if (WeightedSumDenom != 0.0f)
#				{
#					background = WeightedSumBg / WeightedSumDenom;
#					noise = WeightedSumNoise / WeightedSumDenom;
#				}
#
#				BG[iUnit][iAtom].value1 = background;
#				Noise[iUnit][iAtom].value1 = noise;
#
#				float smoothFactor;
#				float scaledAdjustedI;
#				float inten;
#				float modifiedI;
#
#				inten = pCell->GetIntensity((int)Cellxy.x,(int)Cellxy.y);
#				modifiedI = ModifyIntensitySlightly(inten);
#				scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
#				PM[iUnit][iAtom] = scaledAdjustedI;
#				probeSetName = m_Cdf.GetProbeSetName(iUnit);
#				fprintf(ff,"Probeset\t%s\tiUnit\t%i\tiAtom\t%i\tPM\tinten\t%16.8f\tmodified\t%16.8f\tscaled\t%16.8f\n",probeSetName.c_str(),iUnit,iAtom,inten,modifiedI,scaledAdjustedI);
#				
#				//////////////////////////////////////////////////////
#				// Compute Mis-match intensities
#				//////////////////////////////////////////////////////
#				Cellxy.x = mmcell.GetX();
#				Cellxy.y = mmcell.GetY();
#				WeightedSumBg = 0.0f;
#				WeightedSumNoise = 0.0f;
#				WeightedSumDenom = 0.0f;
#				background = 0.0f;
#				noise = 0.0f;
#				for (k = 0; k < NumberZones; k++)
#				{
#					WeightedSumBg    += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].background;
#					WeightedSumNoise += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF) * ZonesInfo.pZones[k].noise;
#					WeightedSumDenom += ComputeWeightAtXY(Cellxy.x, Cellxy.y, ZonesInfo.pZones[k].center.x, ZonesInfo.pZones[k].center.y, smoothF); 
#				}
#				if (WeightedSumDenom != 0.0f)
#				{
#					background = WeightedSumBg / WeightedSumDenom;
#					noise = WeightedSumNoise / WeightedSumDenom;
#				}
#
#				BG[iUnit][iAtom].value2 = background;
#				Noise[iUnit][iAtom].value2 = noise;
#
#
#				inten = pCell->GetIntensity((int)Cellxy.x,(int)Cellxy.y);
#				modifiedI = ModifyIntensitySlightly(inten);
#				scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
#				MM[iUnit][iAtom] = scaledAdjustedI;
#
#				atomCount++;
#			}
#			PM[iUnit].resize(atomCount);
#			MM[iUnit].resize(atomCount);
#		} // if expression type 
#	} // for each unit
#
#	//Clean up
#	delete [] NumberCellsPerZone;
# }
#######################################################################
