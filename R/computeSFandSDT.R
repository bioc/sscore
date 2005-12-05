computeSFandSDT <- function(afbatch,TGT=500, digits=NULL,verbose=FALSE,
 plot.histogram=FALSE,celfile.path=NULL) {
#######################################################################
#
# This function computes the scaling factor (SF) and statistical
# difference threshold (SDT) values of Affymetrix GeneChips, for use in
# calculating S-Score values 
#
# Arguments:
#       afbatch - an AffyBatch object containing the data on the chips
#                 to be compared
#	TGT - the target intensity to which the arrays should be scaled
#       digits - the number of significant digits to retain when
#                returning SF and SDT values
#       verbose - logical value indicating whether additional
#                 information on the calculations should be displayed
#       plot.histogram - logical value indicating whether plots of the
#                        histograms should be displayed
#       celfile.path - character string giving the directory path to
#                      the *.CEL files being read, or NULL if the *.CEL
#                      files are in the current directory
#
# Value:
#       a list containing the components SF (scaling factor values, one
#       for each GeneChip) and SDT (statistical difference threshold
#       values, one for each GeneChip)
#
#######################################################################

# set various parameter values
	c <- 5 # m_Params.TuningConstantCSB
	epsilon <- 0.0001 # m_Params.EpsilonSB

	ContrastTau <- 0.03 # m_Params.ContrastTau
	delta <- 2^-20 # m_Params.Delta
	correction <- 1 + 0 # m_Params.BiasCorrect

	NumberVertZones <- 4
	NumberHorZones <- 4
	Epsilon <- 0.5
	NumberBGCells <- 2
	SmoothFactorBG <- 100
	NoiseFrac <- 0.5
	m.PercentBGCells <- 2

# get the information about the AffyBatch object
	probenames <- geneNames(afbatch)
	intens <- intensity(afbatch)
	nProbeSet <- length(probenames)
	numSamples <- length(afbatch)
	pm.index <- pmindex(afbatch)
	mm.index <- mmindex(afbatch)
	probe.index <- c(unlist(pm.index),unlist(mm.index))
	probe.xy <- indices2xy(probe.index,abatch=afbatch)
	probe.list <- names(probe.index)

# calculate the zones for each probe
	CellsRemaining <- ncol(afbatch) %% NumberVertZones
	if (CellsRemaining != 0 ) 
		zonex <- (ncol(afbatch) + NumberVertZones - CellsRemaining) / NumberVertZones else
		zonex <- ncol(afbatch) / NumberVertZones
	
	CellsRemaining <- nrow(afbatch) %% NumberHorZones
	if (CellsRemaining != 0)
		zoney <- (nrow(afbatch) + NumberHorZones - CellsRemaining) / NumberHorZones else
		zoney <- nrow(afbatch) / NumberHorZones
	
	zoney <- zoney + (zoney %% 2)

	NumberZones <- NumberVertZones * NumberHorZones

	NumUnits <- length(geneNames(afbatch))
	totalCells <- length(intens)

	NumberCellsPerZone <- rep(0,NumberZones)

	fZx <- probe.xy[,"x"] / zonex
	fZy <- probe.xy[,"y"] / zoney
	Zx <- floor(fZx)
	Zy <- floor(fZy)
	probe.zoneID <- Zx + Zy * NumberVertZones

	iZone <- 0:(NumberZones-1)
	x1 <- iZone %% NumberVertZones * zonex
	y1 <- iZone %/% NumberVertZones * zoney
	x2 <- x1 + zonex
	y2 <- y1 + zoney
	ZonesInfo.center.x <- (x1 + x2) / 2
	ZonesInfo.center.y <- (y1 + y2) / 2

	numer <- ncol(afbatch) * nrow(afbatch) * m.PercentBGCells
	denom <- NumberVertZones * NumberHorZones * 100
	bgCells <- floor(numer / denom)

# calculate the rawQ for each *.CEL file
	fname <- sampleNames(afbatch)
	if (length(fname) > 1) 
		rawQ <- sapply(1:length(fname),function(x) computeRawQ(fname[x],intens[,x],probe.index,probe.zoneID,bgCells,NumberZones,celfile.path)) else
		rawQ <- computeRawQ(fname,intens,probe.index,probe.zoneID,bgCells,NumberZones,celfile.path)

# compute the zone information (background and noise for each zone)
	ZoneI <- intens[probe.index,]
	ZoneI <- ifelse(ZoneI > Epsilon,ZoneI,Epsilon)
	if (length(fname) > 1) 
		ZonesInfo <- apply(ZoneI,2,function(x) sapply(split(x,probe.zoneID),computeZoneIInfo,NumberBGCells)) else
		ZonesInfo <- sapply(split(ZoneI,probe.zoneID),computeZoneIInfo,NumberBGCells)
	ZonesInfo.background <- matrix(split(ZonesInfo,c(1,2))[[1]],nrow=NumberVertZones*NumberHorZones)
	ZonesInfo.noise <- matrix(split(ZonesInfo,c(1,2))[[2]],nrow=NumberVertZones*NumberHorZones)

	smoothF <- SmoothFactorBG

# compute the weighted background and noise, which accounts for
# the distance between zones
	Cellxy.x <- probe.xy[,"x"]
	Cellxy.y <- probe.xy[,"y"]
	
	len.xy <- length(Cellxy.x)
	len.zones <- length(ZonesInfo.center.x)
	Denom <- matrix(1/((rep(Cellxy.x,each=len.zones)-rep.int(ZonesInfo.center.x,len.xy))^2 + (rep(Cellxy.y,each=len.zones)-rep.int(ZonesInfo.center.y,len.xy))^2 + smoothF),nrow=len.zones,ncol=len.xy,byrow=FALSE)
	WeightedSumBg <- t(ZonesInfo.background) %*% Denom
	WeightedSumNoise <- t(ZonesInfo.noise) %*% Denom
	WeightedSumDenom <- matrix(1,nrow=length(afbatch),ncol=NumberVertZones*NumberHorZones) %*% Denom
	
	background <- ifelse(WeightedSumDenom != 0,WeightedSumBg / WeightedSumDenom,0)
	noise <- ifelse(WeightedSumDenom != 0,WeightedSumNoise / WeightedSumDenom,0)
	
# find the scaled adjusted intensities
	modifiedI <- pmax(ZoneI,Epsilon)
	factoredNoise <- t(noise * NoiseFrac)
	diff <- modifiedI - t(background)
	scaledAdjustedI <- pmax(diff,factoredNoise,0.5)
		
# split the intensities into PM and MM values
	half <- length(Cellxy.x) %/% 2
	PM <- as.vector(scaledAdjustedI[1:half,])
	MM <- as.vector(scaledAdjustedI[(half+1):length(Cellxy.x),])

# set up the factors for tapply
	probe.group <- sapply(pm.index,length)
	probelist.group <- rep(1:(length(probe.group)*numSamples),rep(probe.group,numSamples)) 

# find the parameters as described in the Affymetrix
# statistical algorithms document
	logPM.minus.logMM <- logb(PM,2) - logb(MM,2)
	SB <- tapply(logPM.minus.logMM,probelist.group,OneStepBiweightAlgorithm,c,epsilon)

	CT.denom <- rep(ifelse(SB > ContrastTau,2^SB,2^(ContrastTau / (1 + (ContrastTau - SB) / 10))),rep(probe.group,numSamples))
	CT <- ifelse(MM < PM,MM,PM / CT.denom)

	v <- PM - CT
	PV <- ifelse(v < delta,correction*logb(delta,2),correction*logb(v,2))

	avgMeasurement <- tapply(PV,probelist.group,function(x) 2^OneStepBiweightAlgorithm(x,c,epsilon))
	avgMeasurement <- as.matrix(avgMeasurement)
	dim(avgMeasurement) <- c(nProbeSet,numSamples)

# now get the SF and SDT
	SF <- TGT / apply(avgMeasurement,2,trimMean,0.02,0.98)
	SDT <- 4 * rawQ * SF

# round to the appropriate number of digits, if desired
# by the user
	if (!is.null(digits)) {
		SF <- round(SF,digits)
		SDT <- round(SDT,digits)
	}

# output the verbose information, if desired by the user
	if (verbose) {
		chip <- whatcdf(fname[1])
		info <- rbind(SF,SDT,rawQ,apply(background,1,mean),sqrt(apply(background,1,var)),
		  apply(background,1,min),apply(background,1,max),apply(noise,1,mean),sqrt(apply(noise,1,var)),
		  apply(noise,1,min),apply(noise,1,max))
		colnames(info) <- fname
		rownames(info) <- c("SF","SDT","rawQ","Background Avg","Background Stdev","Background Min","Background Max",
		  "Noise Avg","Noise Stdev","Noise Min","Noise Max")
		writeLines("SF and SDT parameters:") 
		writeLines(sprintf("Probearray type:\t%s", chip)) 
		writeLines(sprintf("Number of samples:\t%s",numSamples))
		writeLines(sprintf("Number of probesets:\t%d",nProbeSet))
		writeLines(sprintf("Target scale value:\t%d",TGT))
		writeLines(" ")
		print(info)
		writeLines(" ")
	}

# plot the histograms, if desired by the user
	if (plot.histogram) {
		BINw <- 20.0
		histogram <- rep(0,1000) 
		histoPM <- rep(0,1000) 
		rawMM <- intens[unlist(mm.index),]
		if (any(rawMM < 0))
			stop("Raw MM intensity is out of range.  Please check CEL file")
		logMM <- as.matrix(round(BINw * logb(rawMM,2)))
		numMM <- apply(logI,2,table) 

		rawPM <- intens[unlist(pm.index),]
		if (any(rawPM < 0))
			stop("Raw PM intensity is out of range.  Please check CEL file")
		logPM <- as.matrix(round(BINw * logb(rawPM,2))) 
		numPM <- apply(logI,2,table) 

		for (i in 1:ncol(logPM)) {
			x11()
			par(mfrow=c(1,2))
			hist(logPM[,i],breaks=length(numPM[[i]]),main=paste("Histogram of",fname[i]),
				xlab="log PM intensity")
			hist(logMM[,i],breaks=length(numMM[[i]]),main=paste("Histogram of",fname[i]),
				xlab="log MM intensity")
		}
	}

	return(list(SF=SF,SDT=SDT))
}

#######################################################################
# The corresponding Affymetric C++ code follows below
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
#				blk.GetCell(iCell, cell);
#				bMasked = pCell->IsMasked(cell.GetX(), cell.GetY());
#				if (bMasked == false)
#				{
#					int nZone;
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
#				}
#				if (WeightedSumDenom != 0.0f)
#				{
#					background = WeightedSumBg / WeightedSumDenom;
#					noise = WeightedSumNoise / WeightedSumDenom;
#				}
#
#				BG[iUnit][iAtom].value1 = background;
#				Noise[iUnit][iAtom].value1 = noise;
#
#				//float smoothFactor;
#				float scaledAdjustedI;
#				float inten;
#				float modifiedI;
#
#				inten = pCell->GetIntensity((int)Cellxy.x,(int)Cellxy.y);
#				modifiedI = ModifyIntensitySlightly(inten);
#				scaledAdjustedI = ComputeAdjustedIntensity(modifiedI, background, noise);
#				PM[iUnit][iAtom] = scaledAdjustedI;
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
# float computeSquaredDistance(float x1, float y1, float x2, float y2)
# {
# 	float diffx = x1 - x2;
# 	float diffy = y1 - y2;
# 	return (float) ( diffx * diffx + diffy * diffy );
# }
#######################################################################
# float ComputeWeightAtXY(float x, float y, float centerX, float centerY, float smoothFactor)
# {
# 	return 1.0 / (computeSquaredDistance(x , y, centerX, centerY) + smoothFactor);
# }
#######################################################################
# void CExpressionAlgorithmImplementation::ComputeMeasurement(vector<vector<float> > & PM,
# 								  vector<vector<float> > & MM,
# 								  vector<vector<bool> > & UseAtom,
# 								  vector<float> & avgMeasurement,
# 								  vector<vector<float> > & PV)
# {
# 	// Compute the typical difference (for each probe set) of log(PM) over log(MM).
# 	int nProbeSet = PM.size(); //  Number of Probe Set
# 
# 	// Compute Contrast Value
# 	vector<vector<float> > CT(nProbeSet);
# 	ComputeContrastValue(PM, MM, CT);
# 
# 	// Compute Probe Value
# 	ComputeProbeValue(PM, CT, PV);
# 
# 	// Compute Average Log Intensity for probe set i and its confidence interval.
# 	float c = m_Params.TuningConstantCAvgLogInten;
# 	float epsilon = m_Params.EpsilonAvgLogInten;
# 
# 	vector<vector<float> > PVused(nProbeSet);
# 	GetUsedSet(PV, UseAtom, PVused);
# 
# 	vector<float> uncertainty(nProbeSet);
# 	for (int i=0; i<nProbeSet; i++)
# 	{
# 		avgMeasurement[i] = OneStepBiweightAlgorithm(PVused[i], c, epsilon);
# 		avgMeasurement[i] = antiLog(avgMeasurement[i]);
# 	}
# }
#######################################################################
# void CExpressionAlgorithmImplementation::ComputeContrastValue(vector<vector<float> > & PM,
# 									vector<vector<float> > & MM,
# 									vector<vector<float> > & CT)
# {
# 	int nProbeSet = PM.size();
# 	vector<float> SB(nProbeSet);
# 	ComputeTypicalDifference(PM, MM, SB);
# 	float ContrastTau = m_Params.ContrastTau;
# 	for (int i=0; i<nProbeSet; i++)
# 	{
# 		int nProbePair = PM[i].size();
# 		CT[i].resize(nProbePair);
# 		for (int j=0; j<nProbePair; j++)
# 		{
# 			if (MM[i][j] < PM[i][j])
# 			{
# 				CT[i][j] = MM[i][j];
# 			}
# 			else if ((MM[i][j] >= PM[i][j]) &&
# 					 (SB[i] > ContrastTau))
# 			{
# 				CT[i][j] = PM[i][j] / antiLog(SB[i]);
# 			}
# 			else if ((MM[i][j] >= PM[i][j]) &&
# 					 (SB[i] <= ContrastTau))
# 			{
# 				CT[i][j] = PM[i][j] / antiLog(ContrastTau / (1.0 + (ContrastTau - SB[i]) / m_Params.ScaleTau) );
# 			}
# 		}
# 	}
# }
#######################################################################
# void CExpressionAlgorithmImplementation::ComputeProbeValue(vector<vector<float> > & PM,
# 								 vector<vector<float> > & CT,
# 								 vector<vector<float> > & PV)
# {
# 	int nProbeSet = PM.size();
# 	float delta = m_Params.Delta;
# 	float correction = 1.0f + m_Params.BiasCorrect;
# 	for (int i=0; i<nProbeSet; i++)
# 	{
# 		int nProbePair = PM[i].size();
# 		PV[i].resize(nProbePair);
# 		for (int j=0; j<nProbePair; j++)
# 		{
# 			float v = PM[i][j] - CT[i][j];
# 			if (v < delta)
# 				PV[i][j] = correction * logtwo(delta);
# 			else
# 				PV[i][j] = correction * logtwo(v);
# 		}
# 	}
# }
#######################################################################
# void CExpressionAlgorithmImplementation::ComputeTypicalDifference(vector<vector<float> > & PM,
#										vector<vector<float> > & MM,
#										vector<float> & SB)
# {
#	int nProbeSet = PM.size();
#	vector<vector<float> > logPM_minus_logMM(nProbeSet);
#	float c = m_Params.TuningConstantCSB;
#	float epsilon = m_Params.EpsilonSB;
#
#	for (int i=0; i<nProbeSet; i++)
#	{
#		int nProbePair = PM[i].size(); // Number of Probe Pair in Probe Set i
#		logPM_minus_logMM[i].resize(nProbePair);
#		for (int j=0; j<nProbePair; j++)
#		{
#			logPM_minus_logMM[i][j] = logtwo(PM[i][j]) - logtwo(MM[i][j]);
#		}
#		SB[i] = OneStepBiweightAlgorithm(logPM_minus_logMM[i], c, epsilon);
#	}
# }
#######################################################################
# float CExpressionAlgorithmImplementation::DetermineScaleFactor(vector<AbsStatExpressionProbeSetResultType *> &statResults)
# {
#	int iUnit;
#	float avg = 0.0f;
#
#	
#	// User defined norm factor is already stored in the algorithm parameters.
#	if (m_Params.SFMethod == CExpStatAlgSettings::DEFINED_SCALING_FACTOR)
#		return m_Params.ScaleFactor;
# 
#	// Loop over all of the units.
#	int UnitsPerChip = m_Cdf.GetHeader().GetNumProbeSets();
#	int unitCount=0;
#	vector<float> intensityList(UnitsPerChip);
#	string probeSetName;
#	for (iUnit=0; iUnit<UnitsPerChip; iUnit++)
#	{
#		// Get the units.
#		AbsStatExpressionProbeSetResultType *pUnitResult = statResults[iUnit];
#
#		// find signal only for the used genes.
#		probeSetName = m_Cdf.GetProbeSetName(iUnit);
#		if (m_Params.SFMethod == CExpStatAlgSettings::SCALE_TO_ALL_PROBE_SETS ||
#			UseUnitInNormFactor(probeSetName, m_Params.ScaleGenes))
#		{
#			// Use Measurement to do scale factor
#			if (pUnitResult->NumUsedPairs != 0)
#			{
#				intensityList[unitCount] = pUnitResult->Signal;
#				unitCount++;
#			}
#		}
#	}
#
#	intensityList.resize(unitCount);
#
#	// Compute the trimMean
#	float p1 = m_Params.IntensityLowPercent / 100;
#	float p2 = 1.0f - m_Params.IntensityHighPercent / 100;
#	avg = trimMean(intensityList, p1, p2);
#
#	// Store the scale factor
#	float sf=1.0f;
#	if (unitCount && avg != 0.0f)
#		sf = m_Params.TGT / avg;
#
#	//Check for the validity of SF
#	if(sf <= 0)
#	{
#		sf = 1.0f;
#	}
#	return sf;
# }
#######################################################################
# float CExpressionAlgorithmImplementation::ModifyIntensitySlightly(float intensity)
# {
#	return max(intensity, m_Params.Epsilon);
# }
#######################################################################
# float CExpressionAlgorithmImplementation::ComputeAdjustedIntensity(float intensity, float background, float noise)
# {
#	float factoredNoise = (float) noise * m_Params.NoiseFrac;
#	float diff = intensity - background;
#	float adjustedI = max(max(diff, factoredNoise), 0.5f);
#	// AlexC - 1/22/01
#	// Code Comments:
#	// if too frequent substitution of the noise value, alert might generated.
#	// Refer to the page 4 of the Background and Spatial Variation Adjustement spec., 
#	// under eq(16), it said that 
#	// "Production software should probably alert the user if the noise value is being 
#	//  substituted too frequently (indicating that too much data is below the noise level), 
#	//  but an appropriate threshold value is not at present known."
#	return adjustedI;
# }
#######################################################################
# int CExpressionAlgorithmImplementation::DetermineZone(int cellx, int celly, int zonex, int zoney)
# {
#	float fZx = (float) cellx / (float) zonex;
#	float fZy = (float) celly / (float) zoney;
#
#	int Zx = floor(fZx);
#	int Zy = floor(fZy);
#
#	int zoneID = Zx + Zy * (int) m_Params.NumberVertZones;
#	return zoneID;
# }
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
