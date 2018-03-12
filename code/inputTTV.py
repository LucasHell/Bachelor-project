import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pathlib2

G = 0.000295994511
inputFile = []
data = []
name = [] 		# KOI name of planet
kepID = [] 		# kepler ID of planet
mStar = [] 		# mass of host star
mPlanet = [] 	# mass of planet
period = []		# period of planet
rPlanet = []	# radius of planet
numEpoch = []	# number of epochs
eccentricity = 0	# eccentricity of planet
inclination = 90	# icnlination of planet
lNode = 0			# long node of planet
argument = 0		# argument of planet
meanAnom = []		# mean anomaly of planet
count = 0			# counter for number of planets
kepMag = []			# Kepler magnitude of star
transitDur = []		# Duration of transit
rStar = []			# Radius of star
RA = []				# Right ascension of system
dec = []			# Declination of system





#~ with open('nasaDATA.csv','r') as inputFile: # read in data from csv file to respective arrays
	#~ data = inputFile.readlines()[46:]
	#~ name = [k.split(',')[2] for k in data]
	#~ kepID = [k.split(',')[1] for k in data]
	#~ period = [k.split(',')[5] for k in data]
	#~ rPlanet = [k.split(',')[20] for k in data]
	#~ mStar = [k.split(',')[32] for k in data]
	#~ numEpoch = [k.split(',')[8] for k in data]
	#~ kepMag = [k.split(',')[37] for k in data]
	#~ transitDur = [k.split(',')[14] for k in data]
	#~ rStar = [k.split(',')[29] for k in data]
	
	
# read in data from csv file to respective arrays
with open('TESSData.csv','r') as inputFile: 
	data = inputFile.readlines()
	name = [k.split(' ')[0] for k in data]
	kepID = [k.split(' ')[1] for k in data]
	period = [k.split(' ')[2] for k in data]
	rPlanet = [k.split(' ')[3] for k in data]
	mStar = [k.split(' ')[4] for k in data]
	numEpoch = [k.split(' ')[5] for k in data]
	transitDur = [k.split(' ')[6] for k in data]
	rStar = [k.split(' ')[7] for k in data]
	RA = [k.split(' ')[8] for k in data]
	dec = [k.split(' ')[9] for k in data]
	effTemp = [k.split(' ')[10] for k in data]
	kepMag = [k.split(' ')[12] for k in data]
	

#Convert planet radius to mass
for i in range(len(rPlanet)):
	if not rPlanet[i]:
		rPlanet[i] = 0
	else:
		rPlanet[i] = float(rPlanet[i]) 										# convert strings in rPlanet to floats
	
	if rPlanet[i] < 1.5: 													# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(0.440*(rPlanet[i]**3) + 0.614*(rPlanet[i]**4))
	else:																	# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(2.69*(rPlanet[i]**0.93))

	mPlanet[i] = mPlanet[i]*0.000002988 									#convert planet mass from earth masses to solar masses

	
pathlib2.Path('./input/').mkdir(parents=True, exist_ok=True)				# create directory input if it does not exist


# open output files to write data for TTVFast
systemCount = 0	
outputFile = open('input/%s.in' % systemCount, 'w')
outNumP = open('numberPlanets.csv', 'w')
outErrorFile = open('timingErrors.csv', 'w')
outTESSTime = open('TESSTime.csv', 'w')
#~ outEcOm = open('eccenOmeg.csv', 'w')
outputFile.write(repr(G) + '\n' + mStar[0] + '\n')

# write to output files in the format required for TTVFast
for i in range(len(mPlanet)): 
	meanAnom.append(90 - 360 * (float(numEpoch[i]) / float(period[i])))		# calculate mean anomaly of planet with reference to transit of first planet
	while meanAnom[i] > 360:
		meanAnom[i] = meanAnom[i] - 360 									# if angle of mean anomaly is above 360 degrees, subtract 360 until it is in the range 0 to 360 degrees
		if meanAnom[i] == 360:
			meanAnom[i] = 0
				
	while meanAnom[i] < 0:													# if angle of mean anomaly is below 360 degrees, add 360 until it is in the range 0 to 360
		meanAnom[i] = meanAnom[i] + 360
		if meanAnom[i] == 360:
			meanAnom[i] = 0
		
	# write data for TTVFast
	outputFile.write(repr(mPlanet[i]) + '\n' + period[i] + ' ' + str(np.random.rayleigh(0.03)) + ' ' + repr(inclination) + ' ' + repr(lNode) + '  ' + str(np.random.uniform(0,360)) + ' ' + repr(meanAnom[i]) + '\n') 
	count += 1			# counter for number of planets for each system
	
	if i == len(mPlanet)-1: 			# for last element
		if count == 1:
			os.remove('input/%s.in' % systemCount)
			outputFile.close()
			
		elif count != 1:
			outNumP.write(repr(count) + '\n')
			#~ outEcOm.write(str(np.random.rayleigh(0.03)) + ',' + str(np.random.uniform(0,360))+ '\n')
			S = 3.96 * 10**13 * 10**(-0.4*(float(kepMag[i])))
			errorTiming = ((S * float(transitDur[i]))**(-1/2) * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**(-3/2) * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60)) 		# write errors to file in minutes
			outTESSTime.write(RA[i] + ',' + dec[i])	

	elif kepID[i] != kepID[i+1]:		# if ID is not the same as ID of next planet a new file is created for new system
		if count == 1:					# if number of planets is one the file is removed
			os.remove('input/%s.in' % systemCount)
			outputFile.close()
			
			# begin on new file
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			count = 0
			
		elif count != 1:				# if number of planets is not 1 the rest of the data is saved and a new file is created
			outNumP.write(repr(count) + '\n')
			count = 0
			#~ outEcOm.write(str(np.random.rayleigh(0.03)) + ',' + str(np.random.uniform(0,360)) + '\n')
			S = 3.96 * 10**13 * 10**(-0.4*(float(kepMag[i])))
			errorTiming = ((S * float(transitDur[i]))**(-1/2) * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**(-3/2) * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			outTESSTime.write(RA[i] + ',' + dec[i] + '\n')		
			outputFile.close()
			
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			
			
	elif kepID[i] == kepID[i+1] and int(name[i][-2:]) > int(name[i+1][-2:]):
		if count == 1:					# if number of planets is one the file is removed
			os.remove('input/%s.in' % systemCount)
			outputFile.close()
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			count = 0
			
		elif count != 1:				# if number of planets is not 1 the rest of the data is saved and a new file is created
			outNumP.write(repr(count) + '\n')
			count = 0
			#~ outEcOm.write(str(np.random.rayleigh(0.03)) + ',' + str(np.random.uniform(0,360)) + '\n')
			S = 3.96 * 10**13 * 10**(-0.4*(float(kepMag[i])))
			errorTiming = ((S * float(transitDur[i]))**(-1/2) * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**(-3/2) * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			outTESSTime.write(RA[i] + ',' + dec[i] + '\n')
			outputFile.close()
			
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			
		
		
		
		

outNumP.close()
outErrorFile.close()
outTESSTime.close()


	
	


	




