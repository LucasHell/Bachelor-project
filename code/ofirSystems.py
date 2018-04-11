import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


KOINum = []
detFlag = []
TTVFreq = []
ofirP = []
ofirName = []
ofirPlanetR = []
ofirMStar = []
ofirEpochNum = []
ofirMag = []
ofirTransitDur = []
ofirStarR = []



data = pd.read_table('ofir_table.txt', sep=';', skiprows=34, names=('KOI_num', 'newDetFlag', 'TTVfre', 'TTV+uncer', 'TTV-uncer', 'TTV_per', 'Delta_chi', 'chi_area', 'chi_single', 'chi_RMS', 'cho_correl', 'TTV_amp', 'TTV_amp+_uncer', 'TTV_amp-_uncer', 'TTV_ref', 'TTV_ref+_uncer', 'TTV_ref-_uncer', 'cofid', 'STD_error', '20', '21', '22', '23', '24'))
KOINum = data['KOI_num']
detFlag = data['newDetFlag']
TTVFreq = data['TTVfre']
amp_ofir = data['TTV_amp']

with open('nasaDATA.csv','r') as inputFile:
	data = inputFile.readlines()[46:]
	name = [k.split(',')[2] for k in data]
	kepID = [k.split(',')[1] for k in data]
	period = [k.split(',')[5] for k in data]
	rPlanet = [k.split(',')[20] for k in data]
	mStar = [k.split(',')[32] for k in data]
	numEpoch = [k.split(',')[8] for k in data]
	kepMag = [k.split(',')[37] for k in data]
	transitDur = [k.split(',')[14] for k in data]
	rStar = [k.split(',')[29] for k in data]
	
	

outputFileOfir = open('ofirData.csv', 'w')

for i in range(len(KOINum)):
	lenKOI = len(str(KOINum[i]))
	#~ print name[0][lenKOI:]
	for k in range(len(name)):
		if str(KOINum[i]) == str(name[k][-lenKOI:]):
			#~ print KOINum[i], '=', name[k][-lenKOI:]
			outputFileOfir.write(name[k] + ',' +  kepID[k] + ',' + rPlanet[k] + ',' + period[k] + ',' + mStar[k] + ','  + numEpoch[k] + ','  + transitDur[k] + ',' + rStar[k] + ',' + kepMag[k])
			#~ ofirP.append(period[k])
			#~ ofirName.append(name[k])
			#~ ofirPlanetR.append(rPlanet[k])
			#~ ofirMStar.append(mStar[k])
			#~ ofirEpochNum.append(numEpoch[k])
			#~ ofirMag.append(kepMag[k])
			#~ ofirTransitDur.append(transitDur[k])
			#~ ofirStarR.append(rStar[k])
			
			
			
			
			
			
			
			
			
			
