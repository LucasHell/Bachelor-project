import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pathlib2

G = 0.000295994511
inputFile = []
data = []
name = [] 		# kepler name of planet
kepID = [] 		# kepler ID of planet
mStar = [] 		# mass of host star
mPlanet = [] 	# mass of planet
period = []		# period of planet
rPlanet = []	# radius of planet
numPlanets = []	# number of planets
eccentricity = 0	# eccentricity of planet
inclination = 90	# icnlination of planet
lNode = 0			# long node of planet
argument = 0		# argument of planet
meanAnom = []		# mean anomaly of planet
count = 0			# counter



with open('q1_q17_dr25_koi_CONFIRMED.csv','r') as inputFile: # read in data from csv file to respective arrays
	data = inputFile.readlines()[46:]
	name = [k.split(',')[2] for k in data]
	kepID = [k.split(',')[1] for k in data]
	period = [k.split(',')[7] for k in data]
	rPlanet = [k.split(',')[22] for k in data]
	mStar = [k.split(',')[27] for k in data]
	numPlanets = [k.split(',')[28] for k in data]

numPlanets = list(map(int,numPlanets)) 		# convert strings in numPlanets array to integers

for i in range(len(rPlanet)):
	rPlanet[i] = float(rPlanet[i]) 											# convert strings in rPlanet to floats
	if rPlanet[i] < 1.5: 													# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(0.440*(rPlanet[i]**3) + 0.614*(rPlanet[i]**4))
	else:																	# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(2.69*(rPlanet[i]**0.93))

	mPlanet[i] = mPlanet[i]*0.000002988 									#convert planet mass from earth masses to solar masses

	
pathlib2.Path('./input').mkdir(parents=True, exist_ok=True)
	
outputFile = open('input/%s.in' % name[0], 'w')
outputFile.write(repr(G) + '\n' + mStar[0] + '\n')
periodRef = float(period[0])		# period of first planet in system, used as reference for calculating mean anomaly
countLimit = numPlanets[0]
for i in range(len(mPlanet)): 
	if count != 0:
		meanAnom.append(360 * (float(periodRef)/float(period[i])) + 90)		# calculate mean anomaly of planet with reference to transit of first planet
		while meanAnom[i] > 360:
			meanAnom[i] = meanAnom[i] - 360 							# if angle of mean anomaly is above 360 degrees, subtract 360 to keep the range in 0 to 360 degrees
			if meanAnom[i] == 360:
				meanAnom[i] = 0

	else:
		meanAnom.append(90)
		
			
	outputFile.write(repr(mPlanet[i]) + '\n' + period[i] + ' ' + repr(eccentricity) + ' ' + repr(inclination) + ' ' + repr(lNode) + '  ' + repr(argument) + ' ' + repr(meanAnom[i]) + '\n') # + meanAnom
	#~ print count, meanAnom[i]
	count += 1
	
	if name[i] == name[-1]:
		outputFile.close()
	elif kepID[i] != kepID[i+1]:
		if count == 1:
			os.remove('input/%s.in' % name[i])
		outputFile.close()
		outputFile = open('input/%s.in' % name[i+1], 'w')
		outputFile.write(repr(G) + '\n' + mStar[0] + '\n')
		count = 0
		periodRef = period[i+1]
        

#~ for x in range(0,3):
	#~ print name[x], period[x], rPlanet[x], mPlanet[x],  '\n'




