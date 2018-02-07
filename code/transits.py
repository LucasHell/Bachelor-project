import matplotlib.pyplot as plt
import numpy as np
import sys

timesFile = []
valueArray = []
transitTime1Float = []
epoch1Float = []
count = 0
period = 1.0917340278625494e+01


timesFile = open(sys.argv[1])
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

print "Amplitude:", transitAmplitude, "minutes or", transitAmplitude/60, "hours"

if transitMax < 0:
	transitTime1Corrected = transitTimesLinFitted + abs(transitCorrection)
else:
	transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)


plt.plot(epoch1Float*fitTimes[0]/1440, transitTime1Corrected, label='Transit Time')
plt.xlabel('Time [Days]')
plt.ylabel('Transit time [Minutes]')
plt.title('Transit Timing variations')
plt.legend()
plt.savefig('plots/' + sys.argv[1] + '.png')

