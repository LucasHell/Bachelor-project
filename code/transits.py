import matplotlib.pyplot as plt
import numpy as np
import sys

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


for k in range(len(valueArray)):
	if planet[k] == '0':
		epoch1Float.append(float(epoch1[k]))
		transitTime1Float.append(float(transitTime1[k]))

epoch1Float = np.array(epoch1Float)
transitTime1Float = np.array(transitTime1Float)

transitTime1Min = transitTime1Float * 1440
fitTimes = np.polyfit(epoch1Float, transitTime1Min, 1)

transitTimesLinFitted = transitTime1Min-fitTimes[0]*epoch1Float
transitMax = np.amax(transitTimesLinFitted)
transitMin = np.amin(transitTimesLinFitted)
transitAmplitude = (transitMax - transitMin) / 2
transitCorrection = (transitMax + transitMin) / 2

outputFile = open('transAmpl.txt', 'w')
outputFile.write(repr(transitAmplitude))

print "Amplitude:", transitAmplitude, "minutes or", transitAmplitude/60, "hours"

if transitMax < 0:
	transitTime1Corrected = transitTimesLinFitted + abs(transitCorrection)
else:
	transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)

with open('timingErrors.csv','r') as inputFile:
	data = inputFile.readlines()[0:]
	errorTiming = float(data[int(sys.argv[2])])


plt.scatter(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, label='Transit Time')
plt.errorbar(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, yerr = errorTiming, linestyle="None")
plt.xlabel('Time [Days]')
plt.ylabel('Transit time [Minutes]')
plt.title('Transit Timing variations')
plt.legend()
plt.tight_layout()
textstr = 'Amplitude=%.2f\nError=%.2f\n'%(transitAmplitude, errorTiming)
plt.figtext(0.75, 0.5, textstr, fontsize=10)
plt.subplots_adjust(right=0.7)
plt.savefig('plots/' + sys.argv[1] + '.png')

