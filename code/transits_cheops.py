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

with open(sys.argv[1],'r') as timesFile:
	valueArray = timesFile.readlines()
	planet = [k.split(' ')[0] for k in valueArray]
	epoch1 = [k.split(' ')[1] for k in valueArray]
	transitTime1 = [k.split(' ')[2] for k in valueArray]
	

with open('CHEOPS_clones/' + sys.argv[1][10:], 'r') as dataFile:
	periodData = dataFile.readlines()
	
epoch1 = map(int, epoch1)
transitTime1 = map(float, transitTime1)
countPlanetZero = 0

#~ with open('RA_dec_sys.csv','r') as inputFile: 
	#~ data = inputFile.readlines()
	#~ RATess = [k.split(',')[0] for k in data]
	#~ decTess = [k.split(',')[1] for k in data]	
	
#~ with open('rPlanet.csv','r') as inputFile: 
	#~ data = inputFile.readlines()
	#~ rPlanet = [k.split(',')[0] for k in data]
	

with open('CHEOPS_numP_P.csv') as inputFile:
	planetCount = [k.split(' ')[0] for k in inputFile]

#~ for n in range(len(planetCount)):
	#~ planetCount[n] = int(planetCount[n][:1])
with open('error_cheops.csv','r') as inputFile:
		data = inputFile.readlines()[0:]
		errorTiming = float(data[int(sys.argv[2])])

if os.stat(sys.argv[1]).st_size == 0:
	for i in range(0,int(planetCount[int(sys.argv[2])])):
		outputFile = open('transAmplCheops.csv', 'a')
		outputFile.write(str(0) + ',' + sys.argv[1][6:-3] + '\n')
		outputFile.close()
		#~ outErrorFile = open('ampErrorCheops.csv', 'a')
		#~ outErrorFile.write(repr(errorTiming) + '\n')
		#~ outErrorFile.close()
		#~ outRAdec = open('RA_dec_p.csv', 'a')
		#~ outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		#~ outRAdec.close()
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
		#~ outErrorFile.close()
		#~ outRAdec = open('RA_dec_p.csv', 'a')
		#~ outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		#~ outRAdec.close()
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
		#~ print epochAmp, transitAmpEpoch
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

	#~ if transitMax < 0:
		#~ transitTime1Corrected = transitTimesLinFitted + abs(transitCorrection)
	#~ else:
	transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)

	
		
	outErrorFile = open('ampErrorCheops.csv', 'a')
	outErrorFile.write(repr(errorTiming) + '\n')
	outErrorFile.close()
	#~ print RATess[int(sys.argv[2])]*-1
	#~ outRAdec = open('RA_dec_p.csv', 'a')
	#~ if float(decTess[int(sys.argv[2])]) > 0:
		#~ outRAdec.write(str(float(RATess[int(sys.argv[2])])) + ',' + str(float(decTess[int(sys.argv[2])][:-1])) + '\n')
	#~ else:
		#~ outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
	#~ outRAdec.close()
	
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
	
	#~ plt.scatter(epochAmp[1:]*fitTimes[0]/1440, transAmplProg, label='Transit Time')
	#~ plt.xlabel('Time [days]')
	#~ plt.ylabel('Amplitude [min]')
	#~ plt.subplots_adjust(right=0.9)
	#~ plt.savefig('plots/ampPlots/' + sys.argv[1] + '_' + str(i) + '.pdf')
	#~ plt.clf()

	
	planetPeriod.append(float(periodData[3 + 2*i][:9]))
	planetAmp.append(transitAmplitude)
	
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

	
	
#~ outPeriodAmp = open('AmplPeriod.csv', 'a')
#~ periodFrac = 0
#~ for l in range(len(planetPeriod)):
	#~ if str(planetPeriod[l]) != str(planetPeriod[1]):
		#~ periodFrac = planetPeriod[l]/planetPeriod[l-1]
	#~ if planetAmp[l] < 200:
		#~ outPeriodAmp.write(str(periodFrac) + ',' + repr(planetAmp[l]) + ',' + str(planetCount[l]) + '\n')
#~ outPeriodAmp.close()

#~ outPeriodAmp = open('AmplPeriodDouble.csv', 'a')
#~ periodFrac = 0
#~ if len(planetPeriodDouble) > 1:
	#~ for l in range(len(planetPeriodDouble)):
		#~ periodFrac = planetPeriodDouble[1]/planetPeriodDouble[0]
		#~ if planetAmpDouble[l] < 200:
			#~ outPeriodAmp.write(str(periodFrac) + ',' + repr(planetAmpDouble[l]) + ',' + str(planetCount[l]) + '\n')
#~ outPeriodAmp.close()




