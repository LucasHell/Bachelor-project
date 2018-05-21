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

with open('RA_dec_sys.csv','r') as inputFile: 
	data = inputFile.readlines()
	RATess = [k.split(',')[0] for k in data]
	decTess = [k.split(',')[1] for k in data]	

with open('numberPlanets.csv') as inputFile:
	planetCount = [k.split(' ')[0] for k in inputFile]

for n in range(len(planetCount)):
	planetCount[n] = int(planetCount[n][:1])

if os.stat(sys.argv[1]).st_size == 0:
	for i in range(0,planetCount[int(sys.argv[2])]+1):
		outputFile = open('transAmpl.csv', 'a')
		outputFile.write('0' + ',' + str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		outputFile.close()
		outRAdec = open('RA_dec_p.csv', 'a')
		outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		outRAdec.close()
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
		outputFile = open('transAmpl.csv', 'a')
		outputFile.write('0' + ',' + str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		outputFile.close()
		outRAdec = open('RA_dec_p.csv', 'a')
		outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
		outRAdec.close()
		print "Amplitude:", '0', "minutes or", '0', "hours"
		transitCount = 0
		continue
		
	epoch1Float = np.array(epoch1Float)
	transitTime1Float = np.array(transitTime1Float)

	transitTime1Min = transitTime1Float * 1440
	fitTimes = np.polyfit(epoch1Float, transitTime1Min, 1)

	transitTimesLinFitted = transitTime1Min-fitTimes[0]*epoch1Float
	transitMax = np.amax(transitTimesLinFitted)
	transitMin = np.amin(transitTimesLinFitted)
	transitAmplitude = (transitMax - transitMin) / 2
	transitCorrection = (transitMax + transitMin) / 2

	outputFile = open('transAmpl.csv', 'a')
	outputFile.write(repr(transitAmplitude) + ',' + str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
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
	
	outRAdec = open('RA_dec_p.csv', 'a')
	outRAdec.write(str(RATess[int(sys.argv[2])]) + ',' + str(decTess[int(sys.argv[2])]))
	outRAdec.close()
	
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
	epoch1Float = []
	transitTime1Float = []
	transitCount = 0






