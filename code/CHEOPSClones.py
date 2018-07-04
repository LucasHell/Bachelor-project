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
mPClonesCHE = []
perClonesCHE = []
mStarClonesCHE = []
mPlanetCHE = []
numEpochClonesCHE = []
kepMagClonesCHE = []
transDurClonesCHE = []
rPlanetClonesCHE = []
rStarClonesCHE = []
numPCHE = []
mPC_sys = []
pC_sys = []
eccen_sys = []
arg_sys = []
meanAnom_sys = []
G = 0.000295994511
numP = []
rPlanet_sys = []


with open('TESSData.csv','r') as inputFile: 
	data = inputFile.readlines()
	name = [k.split(',')[0] for k in data]
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



betaVar = 0
count1 = 0
j = 0
#~ print len(RA), len(rPlanet)	
for i in range(len(RA)): #
	RA[i] = math.radians(float(RA[i]))
	dec[i] = math.radians(float(dec[i]))

	
	betaVar = math.asin(math.cos(tilt)*math.sin(dec[i]) - math.sin(RA[i])*math.cos(dec[i])*math.sin(tilt))
	lamb1 = (math.sin(tilt)*math.sin(dec[i]) + math.sin(RA[i])*math.cos(dec[i])*math.cos(tilt))/math.cos(betaVar)
	lamb2 = (math.cos(RA[i])*math.cos(dec[i]))/math.cos(betaVar)
	lamb.append(math.degrees(math.atan2(lamb2, lamb1)))
	beta.append(math.degrees(betaVar))
	count1 += 1
	
	mPClones.append(mPlanet[i])
	perClones.append(period[i])
	mStarClones.append(mStar[i])
	numEpochClones.append(numEpoch[i])
	kepMagClones.append(kepMag[i])
	transDurClones.append(transitDur[i])
	rPlanetClones.append(rPlanet[i])
	rStarClones.append(rStar[i])
	
	if i == 0 and kepID[i] != kepID[i+1] and int(str(name[i])[-1:]) >= int(str(name[i+1])[-1:]):	# RA[i] != RA[i+1] and dec[i] != dec[i+1]
		#~ print 'i =', i, 'count1 = ', count1
		numP.append(count1)
		count1 = 0
	if i == len(RA)-1:
		#~ print 'i =', i, 'count1 = ', count1
		numP.append(count1)
		break
	if i > 1 and kepID[i] != kepID[i+1] and int(str(name[i])[-1:]) >= int(str(name[i+1])[-1:]):	#RA[i] != RA[i-1] and dec[i] != dec[i-1]
		#~ print 'i =', i, 'count1 = ', count1
		numP.append(count1)
		count1 = 0
	elif i > 1 and kepID[i] == kepID[i+1] and int(str(name[i])[-1:]) >= int(str(name[i+1])[-1:]):	#RA[i] != RA[i-1] and dec[i] != dec[i-1]
		#~ print 'i =', i, 'count1 = ', count1
		numP.append(count1)
		count1 = 0

	#~ if -40 < math.degrees(betaVar) < 40:
		#~ k = i
		
		#~ print i, lamb[-1], math.degrees(betaVar)
		#~ mPClones.append(mPlanet[k])
		#~ perClones.append(period[k])
		#~ mStarClones.append(mStar[k])
		#~ numEpochClones.append(numEpoch[k])
		#~ kepMagClones.append(kepMag[k])
		#~ transDurClones.append(transitDur[k])
		#~ rPlanetClones.append(rPlanet[k])
		#~ rStarClones.append(rStar[k])
		#~ numbP = 0
		#~ count1 += 1
		#~ print i	
		#~ numbP += 1
		#~ k += 1
		#~ if kepID[k] == kepID[-1]:
				#~ break
		#~ while kepID[k] == kepID[k+1]:
			
			#~ if kepID[k] == kepID[-1]:
				#~ break
			#~ numbP += 1
			#~ k += 1
			
				
		#~ numP.append(numbP+1)
		#~ if float(dec[i]) > 0:
			#~ outRADec.write(str(math.degrees(RA[i])) + ',' + str(math.degrees(dec[i])*-1) + '\n')
		#~ else:
			#~ outRADec.write(str(math.degrees(RA[i])) + ',' + str(math.degrees(dec[i])) + '\n')
		#~ if numbP+1 > 1:
			#~ outNumP.write(str(numbP+1) + '\n')
			

		

outNumP = open('CHEOPS_numP.csv', 'w')
outErrorTiming = open('error_cheops.csv', 'w')
outNumPP = open('CHEOPS_numP_P.csv', 'w')
systemCount = 0	

k = 0
for l in range(len(numP)):
	if -40 < beta[k] < 40:
		for o in range(0,numP[l]):
			mPClonesCHE.append(mPlanet[k+o])
			perClonesCHE.append(period[k+o])
			mStarClonesCHE.append(mStar[k+o])
			numEpochClonesCHE.append(numEpoch[k+o])
			kepMagClonesCHE.append(kepMag[k+o])
			transDurClonesCHE.append(transitDur[k+o])
			rPlanetClonesCHE.append(rPlanet[k+o])
			rStarClonesCHE.append(rStar[k+o])
		numPCHE.append(numP[l])
	k += numP[l]
	

l = 0
u = 0
count = 0
countP = 0

for l in range(len(numPCHE)):		#
	if numPCHE[l] == 1:
		#~ print count
		count += 1
		u += 1
		#~ print "skip"
	if numPCHE[l] != 1:
		S = 3.958 * 10**11 * 10**(-0.4*float(kepMagClonesCHE[l]))			
		errorTiming = ((S * float(transDurClonesCHE[l]))**Fraction('-1/2') * ((float(rPlanetClonesCHE[l])*0.009158)/float(rStarClonesCHE[l]))**Fraction('-3/2') * float(transDurClonesCHE[l]))
		#~ print numPCHE[l]
		for k in range(0,100):
			outputFile = open('CHEOPS_clones/%s' % systemCount + '_' + str(k) + '.in', 'w')
			outRPlanet = open('radius_cheops/%s' % systemCount + '_' + str(k) + '.txt', 'w')
			outputFile.write(repr(G) + '\n' + mStar[l] + '\n')
			
			for j in range(0,numPCHE[l]):
				argument = np.random.uniform(0,360)
				meanAnom = 90 - (360 * (float(numEpochClonesCHE[l]) / float(perClonesCHE[l]))) - argument
				while meanAnom > 360:
					meanAnom = meanAnom - 360 									
					if meanAnom == 360:
						meanAnom = 0
					
				while meanAnom < 0:		
					meanAnom = meanAnom + 360
					if meanAnom == 360:
						meanAnom = 0
				#~ print u+j, l
				mPC_sys.append(mPClonesCHE[u+j])
				pC_sys.append(perClonesCHE[u+j])
				eccen_sys.append(np.random.rayleigh(0.03))
				arg_sys.append(argument)
				meanAnom_sys.append(meanAnom)
				rPlanet_sys.append(rPlanetClonesCHE[u+j])
				#~ print rPlanetClonesCHE[u+j]

				outNumPP.write(str(numPCHE[l]) + '\n')
				outErrorTiming.write(str(errorTiming) +'\n')
				pC_sys = map(float, pC_sys)
			
			#~ print u
			for k in range(len(mPC_sys)):
				posMin = pC_sys.index(np.amin(pC_sys))

				
				outputFile.write(str(mPC_sys[posMin]) + '\n' + str(pC_sys[posMin]) + ' ' + str(eccen_sys[posMin]) + ' ' + str(90) + ' ' + str(0) + ' ' + str(arg_sys[posMin]) + ' ' + str(meanAnom_sys[posMin]) + '\n')
				outRPlanet.write(str(rPlanet_sys[posMin]) + '\n')
				#~ print rPlanet_sys[posMin], mPC_sys[posMin]
				pC_sys[posMin] = 1000000
			outputFile.close()
			mPC_sys = []
			pC_sys = []
			eccen_sys = []
			arg_sys = []
			meanAnom_sys = []
			rPlanet_sys = []
		systemCount += 1
		u = l + numPCHE[l]
		outNumP.write(str(numPCHE[l]) + '\n')
		countP += numPCHE[l]
	
outNumP.close()		
#~ print countP		
		
		

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
