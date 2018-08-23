import matplotlib.pyplot as plt
import numpy as np
import sys
import os

timesFile = []
valueArray = []
transitTime1Float = []
epoch1Float = []
count = 0
kepMag = []
transitDur = []
rStar = []
planetPeriod = []
planetAmp = []
transitAmpEpoch = []
transAmplProg = []
transTimeProg = []
epochAmp = []
transitTimesLinFittedProg = []
planetPeriodDouble = []
planetAmpDouble = []
rPlanet = []
rPlanetFrac = []

with open(sys.argv[1],'r') as timesFile:
	valueArray = timesFile.readlines()
	planet = [k.split(' ')[0] for k in valueArray]
	epoch1 = [k.split(' ')[1] for k in valueArray]
	transitTime1 = [k.split(' ')[2] for k in valueArray]
	

with open('CHEOPS_clones/' + sys.argv[1][10:], 'r') as dataFile:
	periodData = dataFile.readlines()

with open('radius_cheops/' + sys.argv[1][10:-2] + 'txt', 'r') as inputFile: 
	rPlanet = inputFile.readlines()
	
epoch1 = map(int, epoch1)
transitTime1 = map(float, transitTime1)
countPlanetZero = 0


	

with open('CHEOPS_numP_P.csv') as inputFile:
	planetCount = [k.split(' ')[0] for k in inputFile]


with open('error_cheops.csv','r') as inputFile:
		data = inputFile.readlines()[0:]
		errorTiming = float(data[int(sys.argv[2])])

if os.stat(sys.argv[1]).st_size == 0:
	for i in range(0,int(planetCount[int(sys.argv[2])])):
		outputFile = open('transAmplCheops.csv', 'a')
		outputFile.write(str(0) + ',' + sys.argv[1][6:-3] + '\n')
		outputFile.close()
		print "No value"
	sys.exit(0)
		
transitCount = 0


for i in range(0, int(max(planet))+1): 
	for l in range(len(planet)):
		if planet[l] == str(i):
			epoch1Float.append(int(epoch1[l]))
			transitTime1Float.append(float(transitTime1[l]))
			transitCount += 1

	if transitCount < 2:
		outputFile = open('transAmplCheops.csv', 'a')
		outputFile.write(str(0) + ',' + sys.argv[1][6:-3] + '\n')
		outputFile.close()
		outErrorFile = open('ampErrorCheops.csv', 'a')
		outErrorFile.write(repr(errorTiming) + '\n')
		print "Amplitude:", '0', "minutes or", '0', "hours"
		transitCount = 0
		continue
	

	epoch1Float = np.array(epoch1Float)
	transitTime1Float = np.array(transitTime1Float)
	if epoch1Float[0] == epoch1Float[1]:
		epoch1Float[1:] += 1

	transitTime1Min = transitTime1Float * 1440

	for l in range(len(transitTime1Float)):
		transitAmpEpoch.append(transitTime1Min[l])
		epochAmp.append(epoch1Float[l])
		if len(transitAmpEpoch) > 1:
			fitTimes = np.polyfit(epochAmp, transitAmpEpoch, 1)
			transitTimesLinFittedProg.append(transitAmpEpoch[l]-fitTimes[0]*epoch1Float[l])
			transitMax = np.amax(transitTimesLinFittedProg)
			transitMin = np.amin(transitTimesLinFittedProg)
			transitCorrection = (transitMax + transitMin) / 2	
			transAmplProg.append(transitMax - transitMin)
			
			
	fitTimes = np.polyfit(epoch1Float, transitTime1Min, 1)

	transitTimesLinFitted = transitTime1Min-fitTimes[0]*epoch1Float
	transitMax = np.amax(transitTimesLinFitted)
	transitMin = np.amin(transitTimesLinFitted)
	transitAmplitude = (transitMax - transitMin)
	transitCorrection = (transitMax + transitMin) / 2


	outputFile = open('transAmplCheops.csv', 'a')
	outputFile.write(repr(transitAmplitude) + ',' + sys.argv[1][6:-3] + '\n')
	outputFile.close()

	print "Amplitude:", transitAmplitude, "minutes or", transitAmplitude/60, "hours"

	transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)

	
		
	outErrorFile = open('ampErrorCheops.csv', 'a')
	outErrorFile.write(repr(errorTiming) + '\n')
	outErrorFile.close()
	epochAmp = np.array(epochAmp)

	y_line = np.linspace(0, 0 , len(epoch1Float))
	plt.scatter(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, label='Transit Time')
	plt.errorbar(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, yerr = errorTiming, linestyle="None")
	plt.axhline(y = 0, xmin = 0, xmax = np.amax(epoch1Float*fitTimes[0]/1440), c = 'black')
	plt.xlabel('Time [Days]')
	plt.ylabel('O-C [Minutes]')
	plt.tight_layout()
	textstr = 'Amplitude=%.2f min\nError=%.2f min\n'%(transitAmplitude, errorTiming)
	plt.figtext(0.76, 0.5, textstr, fontsize=10)
	plt.subplots_adjust(right=0.75)
	plt.savefig('plots/' + sys.argv[1] + '_' + str(i) + '.pdf')
	plt.clf()
	

	
	planetPeriod.append(float(periodData[3 + 2*i][:9]))
	planetAmp.append(transitAmplitude)
	rPlanetFrac.append(rPlanet[i][:-2])
	
	if int(max(planet))+1 == 2:
		planetPeriodDouble.append(float(periodData[3 + 2*i][:9]))
		planetAmpDouble.append(transitAmplitude)
	
	epoch1Float = []
	transitTime1Float = []
	transAmplProg = []
	transitAmpEpoch = []
	epochAmp = []
	transitTimesLinFittedProg = []
	transitCount = 0

	

outPeriodAmp = open('AmplPeriod.csv', 'a')
periodFrac = 0
#~ print planetPeriod
for l in range(len(planetPeriod)):
	if str(planetPeriod[l]) == str(planetPeriod[-1]):
		outPeriodAmp.write(str(float(planetPeriod[l])/(float(planetPeriod[l-1]))) + ',' + str(rPlanet[l][:-2]) + ',' + str(planetAmp[l]) + ',' + str(errorTiming) + '\n')
	elif str(planetPeriod[l]) == str(planetPeriod[0]):
		outPeriodAmp.write(str(float(planetPeriod[l+1])/(float(planetPeriod[l]))) + ',' + str(rPlanet[l][:-2]) + ',' + str(planetAmp[l]) + ',' + str(errorTiming) + '\n')
	else:
		outPeriodAmp.write(str(float(planetPeriod[l])/(float(planetPeriod[l-1]))) + ',' + str(rPlanet[l][:-2]) + ',' + str(planetAmp[l]) + ',' + str(errorTiming) + '\n')
		
		outPeriodAmp.write(str(float(planetPeriod[l+1])/(float(planetPeriod[l]))) + ',' + str(rPlanet[l][:-2]) + ',' + str(planetAmp[l]) + ',' + str(errorTiming) + '\n')
	
	
outPeriodAmp.close()





