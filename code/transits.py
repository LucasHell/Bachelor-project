import matplotlib.pyplot as plt
import numpy as np
import sys
import os

timesFile = []
valueArray = []
transitTime1Float = []
epoch1Float = []
count = 0
period = 1.0917340278625494e+01
kepMag = []
transitDur = []
rStar = []

with open(sys.argv[1],'r') as timesFile:
	valueArray = timesFile.readlines()
	planet = [k.split(' ')[0] for k in valueArray]
	epoch1 = [k.split(' ')[1] for k in valueArray]
	transitTime1 = [k.split(' ')[2] for k in valueArray]

epoch1 = map(int, epoch1)
transitTime1 = map(float, transitTime1)
countPlanetZero = 0

	

with open('numberPlanets.csv') as inputFile:
	planetCount = [k.split(' ')[0] for k in inputFile]

for n in range(len(planetCount)):
	planetCount[n] = int(planetCount[n][:1])

if os.stat(sys.argv[1]).st_size == 0:
	for i in range(0,planetCount[int(sys.argv[2])]+1):
		outputFile = open('transAmpl.txt', 'a')
		outputFile.write('0' + '\n')
		outputFile.close()
		print "No value"
	sys.exit(0)
		

transitCount = 0

for i in range(0, int(max(planet))+1):
	#~ print i
	for l in range(len(planet)):
		if planet[l] == str(i):
			epoch1Float.append(int(epoch1[l]))
			transitTime1Float.append(float(transitTime1[l]))
			transitCount += 1
			#~ print "count ", transitCount
	if transitCount < 2:
		outputFile = open('transAmpl.txt', 'a')
		outputFile.write('0' + '\n')
		outputFile.close()
		print "Amplitude:", '0', "minutes or", '0', "hours"
		transitCount = 0
		continue
		
	#~ print "i = ", i
	epoch1Float = np.array(epoch1Float)
	transitTime1Float = np.array(transitTime1Float)

	transitTime1Min = transitTime1Float * 1440
	fitTimes = np.polyfit(epoch1Float, transitTime1Min, 1)

	transitTimesLinFitted = transitTime1Min-fitTimes[0]*epoch1Float
	transitMax = np.amax(transitTimesLinFitted)
	transitMin = np.amin(transitTimesLinFitted)
	transitAmplitude = (transitMax - transitMin) / 2
	transitCorrection = (transitMax + transitMin) / 2

	outputFile = open('transAmpl.txt', 'a')
	outputFile.write(repr(transitAmplitude) + '\n')
	outputFile.close()

	print "Amplitude:", transitAmplitude, "minutes or", transitAmplitude/60, "hours"

	if transitMax < 0:
		transitTime1Corrected = transitTimesLinFitted + abs(transitCorrection)
	else:
		transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)

	with open('timingErrors.csv','r') as inputFile:
		data = inputFile.readlines()[0:]
		errorTiming = float(data[int(sys.argv[2])])
		
	outErrorFile = open('ampError.txt', 'a')
	outErrorFile.write(repr(errorTiming) + '\n')
	outErrorFile.close()
	
	
	y_line = np.linspace(0, 0 , len(epoch1Float))
	#~ print len(epoch1Float), len(transitTime1Corrected)
	plt.scatter(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, label='Transit Time')
	plt.errorbar(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, yerr = errorTiming, linestyle="None")
	plt.axhline(y = 0, xmin = 0, xmax = np.amax(epoch1Float*fitTimes[0]/1440), c = 'black')
	#~ plt.plot(epoch1Float*fitTimes[0]/1440, y_line, c = 'black')
	plt.xlabel('Time [Days]')
	plt.ylabel('Transit time [Minutes]')
	plt.title('Transit Timing variations')
	plt.legend()
	plt.tight_layout()
	textstr = 'Amplitude=%.2f\nError=%.2f\n'%(transitAmplitude, errorTiming)
	plt.figtext(0.75, 0.5, textstr, fontsize=10)
	plt.subplots_adjust(right=0.7)
	plt.savefig('plots/' + sys.argv[1] + '_' + str(i) + '.png')
	plt.clf()
	epoch1Float = []
	transitTime1Float = []
	transitCount = 0






