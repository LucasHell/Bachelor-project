import math
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from fractions import Fraction  

beta = []
lamb1 = 0
lamb2 = 0
lamb = []
tilt = math.radians(23.439281)
mPClones = []
perClones = []
mStarClones = []
mPlanet = []
numEpochClones = []
kepMagClones = []
transDurClones = []
rPlanetClones = []
rStarClones = []
mPC_sys = []
pC_sys = []
eccen_sys = []
arg_sys = []
meanAnom_sys = []
G = 0.000295994511
numP = []


with open('TESSData.csv','r') as inputFile: 
	data = inputFile.readlines()
	#~ name = [k.split(',')[0] for k in data]
	kepID = [k.split(',')[1] for k in data]
	period = [k.split(',')[2] for k in data]
	rPlanet = [k.split(',')[3] for k in data]
	mStar = [k.split(',')[4] for k in data]
	numEpoch = [k.split(',')[5] for k in data]
	transitDur = [k.split(',')[6] for k in data]
	rStar = [k.split(',')[7] for k in data]
	RA = [k.split(',')[8] for k in data]
	dec = [k.split(',')[9] for k in data]
	#~ effTemp = [k.split(',')[10] for k in data]
	kepMag = [k.split(',')[12] for k in data]
	
for i in range(len(rPlanet)):
	if not rPlanet[i]:
		rPlanet[i] = 0
	else:
		rPlanet[i] = float(rPlanet[i]) 										# convert strings in rPlanet to floats
	
	if rPlanet[i] < 1.5: 													# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(0.440*(rPlanet[i]**3) + 0.614*(rPlanet[i]**4))
	else:																	# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(2.69*(rPlanet[i]**0.93))

	mPlanet[i] = mPlanet[i]*0.000002988 		#convert planet mass from earth masses to solar masses	
	
outRADec = open('CHEOPS_RA_Dec.csv', 'w')
outNumP = open('CHEOPS_numP.csv', 'w')
j = 0
betaVar = 0
count1 = 0
print len(RA), len(rPlanet)	
for i in range(len(RA)):
	RA[i] = math.radians(float(RA[i]))
	dec[i] = math.radians(float(dec[i]))

	if RA[i] == RA[j] and dec[i] == dec[j]:
		continue
	betaVar = math.asin(math.cos(tilt)*math.sin(dec[i]) - math.sin(RA[i])*math.cos(dec[i])*math.sin(tilt))
	lamb1 = (math.sin(tilt)*math.sin(dec[i]) + math.sin(RA[i])*math.cos(dec[i])*math.cos(tilt))/math.cos(betaVar)
	lamb2 = (math.cos(RA[i])*math.cos(dec[i]))/math.cos(betaVar)
	lamb.append(math.degrees(math.atan2(lamb2, lamb1)))
	beta.append(math.degrees(betaVar))
	
	j = i
	if -40 < math.degrees(betaVar) < 40:
		#~ print i, lamb[-1], math.degrees(betaVar)
		mPClones.append(mPlanet[i])
		perClones.append(period[i])
		mStarClones.append(mStar[i])
		numEpochClones.append(numEpoch[i])
		kepMagClones.append(kepMag[i])
		transDurClones.append(transitDur[i])
		rPlanetClones.append(rPlanet[i])
		rStarClones.append(rStar[i])
		k = i
		numbP = 0
		count1 += 1
		#~ print i

		while kepID[k] == kepID[k+1]:
			numbP += 1
			k += 1
			if kepID[k] == kepID[-1]:
				break
		numP.append(numbP+1)
		if float(dec[i]) > 0:
			outRADec.write(str(math.degrees(RA[i])) + ',' + str(math.degrees(dec[i])*-1) + '\n')
		else:
			outRADec.write(str(math.degrees(RA[i])) + ',' + str(math.degrees(dec[i])) + '\n')
		outNumP.write(str(numbP+1) + '\n')

		

print len(numP), len(mPClones), count1
outErrorTiming = open('error_cheops.csv', 'w')
outNumPP = open('CHEOPS_numP_P.csv', 'w')
systemCount = 0	
#~ for l in range(len(mPClones)):
l = 0
#~ print mPClones[-1], mPClones[-2]
while l < len(mPClones)-1:
	#~ print l
	if l + numP[l] > len(mPClones):
			break
	if numP[l] != 1:
		S = 3.958 * 10**11 * 10**(-0.4*float(kepMagClones[l]))			
		errorTiming = ((S * float(transDurClones[l]))**Fraction('-1/2') * ((float(rPlanetClones[l])*0.009158)/float(rStarClones[l]))**Fraction('-3/2') * float(transDurClones[l]))
			#~ if l == len(mPClones)-100:
				#~ for i in range(0,100)
				#~ outputFile = open('CHEOPS_clones/%s' % systemCount + '_' + str(k) + '.in', 'w')
				#~ outputFile.write(repr(G) + '\n' + mStar[l] + '\n')
			
		for k in range(0,100):
			outputFile = open('CHEOPS_clones/%s' % systemCount + '_' + str(k) + '.in', 'w')
			outputFile.write(repr(G) + '\n' + mStar[l] + '\n')
			
			for j in range(0,numP[l]):
				argument = np.random.uniform(0,360)
				meanAnom = 90 - (360 * (float(numEpochClones[l]) / float(perClones[l]))) - argument
				while meanAnom > 360:
					meanAnom = meanAnom - 360 									
					if meanAnom == 360:
						meanAnom = 0
					
				while meanAnom < 0:		
					meanAnom = meanAnom + 360
					if meanAnom == 360:
						meanAnom = 0
			
				mPC_sys.append(mPClones[l+j])
				pC_sys.append(perClones[l+j])
				eccen_sys.append(np.random.rayleigh(0.03))
				arg_sys.append(argument)
				meanAnom_sys.append(meanAnom)
				outNumPP.write(str(numP[l]) + '\n')
				outErrorTiming.write(str(errorTiming) +'\n')
				pC_sys = map(float, pC_sys)
				
			for k in range(len(mPC_sys)):
				posMin = pC_sys.index(np.amin(pC_sys))
				#~ print posMin
				outputFile.write(str(mPClones[posMin]) + '\n' + str(perClones[posMin]) + ' ' + str(eccen_sys[posMin]) + ' ' + str(90) + ' ' + str(0) + ' ' + str(arg_sys[posMin]) + ' ' + str(meanAnom_sys[posMin]) + '\n')
				pC_sys[posMin] = 1000000
			outputFile.close()
			mPC_sys = []
			pC_sys = []
			eccen_sys = []
			arg_sys = []
			meanAnom_sys = []
		systemCount += 1

	l += numP[l]
	
		
		
		
		

#~ c = SkyCoord(lon=lamb*u.deg, lat=beta*u.deg, frame='heliocentrictrueecliptic')
#~ ra_rad = c.lon.wrap_at(180 * u.deg).rad			
#~ dec_rad = c.lat.rad


#~ plt.figure(figsize=(8,4.2))
#~ plt.subplot(111, projection="aitoff")
#~ plt.grid(True)
#~ plt.title("Position of observed TESS objects", y=1.08)
#~ plt.axhspan(math.radians(-40), math.radians(40), facecolor='g', alpha=0.1)
#~ plt.scatter(ra_rad, dec_rad, s=7)
#~ plt.show()
