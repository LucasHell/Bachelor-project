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
numPlanets = []	# number of planets
numEpoch = []	# number of epochs
eccentricity = 0	# eccentricity of planet
inclination = 90	# icnlination of planet
lNode = 0			# long node of planet
argument = 0		# argument of planet
meanAnom = []		# mean anomaly of planet
count = 0			# counter for number of planets
kepMag = []			# Kepler magnitude of star
transitDur = []		# Duration of transit
rStar = []			




with open('nasaDATA.csv','r') as inputFile: # read in data from csv file to respective arrays
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
	
	#numPlanets = [k.split(',')[20] for k in data]

numPlanets = list(map(int,numPlanets)) 		# convert strings in numPlanets array to integers

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

	
pathlib2.Path('./input/').mkdir(parents=True, exist_ok=True)
	
outputFile = open('input/%s.in' % name[0], 'w')
outNumP = open('numberPlanets.csv', 'w')
outErrorFile = open('timingErrors.csv', 'w')
outputFile.write(repr(G) + '\n' + mStar[0] + '\n')

for i in range(len(mPlanet)): 
	#~ numEpoch[i] = float(numEpoch[i])							# convert strings to float in epoch array
	#~ period[i] = float
	meanAnom.append(90 - 360 * (float(numEpoch[i]) / float(period[i])))		# calculate mean anomaly of planet with reference to transit of first planet
	while meanAnom[i] > 360:
		meanAnom[i] = meanAnom[i] - 360 							# if angle of mean anomaly is above 360 degrees, subtract 360 to keep the range in 0 to 360 degrees
		if meanAnom[i] == 360:
			meanAnom[i] = 0
				
	while meanAnom[i] < 0:
		meanAnom[i] = meanAnom[i] + 360
		if meanAnom[i] == 360:
			meanAnom[i] = 0
		
			
	outputFile.write(repr(mPlanet[i]) + '\n' + period[i] + ' ' + repr(eccentricity) + ' ' + repr(inclination) + ' ' + repr(lNode) + '  ' + repr(argument) + ' ' + repr(meanAnom[i]) + '\n') # + meanAnom
	#~ print count, meanAnom[i]
	count += 1
	
	if name[i] == name[-1]: 			# for last element
		if count == 1:
			os.remove('input/%s.in' % name[i])
			outputFile.close()
			
		elif count != 1:
			outNumP.write(repr(count) + '\n')
			errorTiming = ((S * float(transitDur[i]))**(-1/2) * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**(-3/2) * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes

	elif kepID[i] != kepID[i+1]:		# for new system
		if count == 1:					# if number of planets is zero the file is removed
			os.remove('input/%s.in' % name[i])
			outputFile.close()
			outputFile = open('input/%s.in' % name[i+1], 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			
		elif count != 1:				# if number of planets is not zero the rest of the data is saved and a new file is created
			outNumP.write(repr(count) + '\n')
			S = 7.8 * 10**8 * 10**(-0.4*(float(kepMag[i])-12))
			#print i, rPlanet[i], rStar[i]
			errorTiming = ((S * float(transitDur[i]))**(-1/2) * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**(-3/2) * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			outputFile.close()
			outputFile = open('input/%s.in' % name[i+1], 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
		
		
		
		count = 0

outNumP.close()
outErrorFile.close()
	


	




