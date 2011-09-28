deltaT		15
stepLimit	100000
stepAvg		10000
resetTime  	1
doCheckpoint	1			
stepCheckpoint	5000
stepPdb		50000
applyThermostat	1
extTemperature	22.5
tauT 		200
applyBarostat	0
keepSquare	1
keepTetragonal	1
extPressure	1
tauP		500
flexBox         1
nebrTabFac	200
rCutLipLip	1.2 
rCutWatWat	0.9
rNebrShell	0.1
stepNebr	10
removeSystemTranslation	1
removeMonolayersTranslation	1
*** MISCELLANY
randSeed	17 
recordSnap	0 
runId		1 
zReplicas	2
*** INITIAL STRUCTURE GENERATION
loadStructure	
region          12.34 12.34 5.05
adjustRegion	 
regionAdjusted	
nSites          12288
initHalfCellWat 48 48 1
nWaters		4608
nSpecies        11
nTypes		7
nDOPEs		512
nLipids		512
reCenterBilayer 0
