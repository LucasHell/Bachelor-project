import matplotlib.pyplot as plt
import numpy as np

timesFile = []
valueArray = []
transitTime1Float = []
epoch1Float = []
count = 0
period = 1.0917340278625494e+01


timesFile = open('Times')
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
transitAmplitude = transitMax - transitMin
transitCorrection = (transitMax + transitMin) / 2

print "Amplitude:", transitAmplitude, "minutes"

if transitMax < 0:
	transitTime1Corrected = transitTimesLinFitted + abs(transitCorrection)
else:
	transitTime1Corrected = transitTimesLinFitted - abs(transitCorrection)


plt.plot(epoch1Float[0:100]*fitTimes[0]/1440, transitTime1Corrected[0:100], label='Transit Time')
plt.xlabel('Times [Days]')
plt.ylabel('Transit time [Minutes]')
plt.title('Transit Timing variations')
plt.legend()
plt.show()

