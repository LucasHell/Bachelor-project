import matplotlib.pyplot as plt
import numpy as np

P = 150
Tobs = 365
A = 48.46
Ameas = []
Acorr = []
t = np.linspace(0,Tobs, 365)
x = []
transitTimesLinFittedProg1 = []
transAmplProg1 = []
transitTimesLinFittedProg = []
transAmplProg = []
epoch1Float = []
transitTime1Float= []
transitCount = 0
transitAmpEpoch = []
epochAmp = []
transitTime1Min = []

with open('./code/times/247.in','r') as timesFile:
	valueArray = timesFile.readlines()
	planet = [k.split(' ')[0] for k in valueArray]
	epoch1 = [k.split(' ')[1] for k in valueArray]
	transitTime1 = [k.split(' ')[2] for k in valueArray]
	
epoch1 = map(int, epoch1)

for i in range(0, Tobs):
	Ameas.append((A/2)*np.sin((2*np.pi*i)/P)-((A/2)*np.sin((2*np.pi*Tobs)/P))*(i/Tobs))
	if i > 200:
		Acorr.append(Acorr[i-1])
	else:
		Acorr.append((Ameas[i]/Tobs)*i)
	x.append(i)

	if len(Ameas) > 1:
		fitTimes = np.polyfit(x, Ameas, 1)
		transitTimesLinFittedProg1.append(Ameas[i]-fitTimes[0]*i)
		transitMax1 = np.amax(transitTimesLinFittedProg1)
		transitMin1 = np.amin(transitTimesLinFittedProg1)
		transAmplProg1.append(transitMax1 - transitMin1)
	
	
	
	
for l in range(len(planet)):
	if planet[l] == str(1):
		epoch1Float.append(int(epoch1[l]))
		transitTime1Float.append(float(transitTime1[l]))
		transitCount += 1

for i in range(len(transitTime1Float)):
	transitTime1Min.append(transitTime1Float[i] * 1440)

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
		epochAmp[l] = epochAmp[l]*fitTimes[0]/1440
		

for o in range(len(epochAmp)):
	epochAmp[o] = epochAmp[o]*fitTimes[0]/1440
	

plt.scatter(epochAmp[1:], transAmplProg, label='Simulated')
plt.plot(x[1:], transAmplProg1, color='black', linestyle='--', label='Analytical')
plt.xlabel('Time [days]')
plt.ylabel('Amplitude [min]')
plt.legend()
plt.show()
plt.savefig('test.pdf')
plt.clf()



#~ xmin = np.amin(1)
#~ xmax = np.amax(365)
#~ x = np.linspace(xmin, xmax, len(transitTimesLinFittedProg))
#~ plt.plot(x, transitTimesLinFittedProg)
#~ plt.show()
