import matplotlib.pyplot as plt
import numpy as np
import math
y = []
x = []
ampProg = []
A = 40
P = 365/2
tobs = 365
ampMax = 0
ampMin = 0
amp = []
Ameas = []

for i in range(0,365):
	Ameas.append((A/2)*np.sin((2*np.pi*i)/P)-((A/2)*np.sin((2*np.pi*tobs)/P))*(i/tobs))
	y.append(40* math.sin(-2*math.radians(i)))
	x.append(i)

	if len(y) > 1:
		linFit = np.polyfit(x, y, 1)
		#~ print y[i], linFit[0]
		y[i] = y[i] - linFit[0]*i
		#~ ampProg.append(np.amax(y)-np.amin(y))
		#~ plt.plot(x,y)
		#~ plt.savefig('./sinText/' + str(i) + '.png')
		#~ plt.clf()
		ampMax = np.amax(Ameas)
		ampMin = np.amin(Ameas)
		
		amp.append((ampMax+ampMin)/2)
	else: 
		ampProg.append(0)
		
		ampMax = np.amax(ampProg)
		ampMax = np.amin(ampProg)
		
		amp.append((ampMax+ampMin)/2)




plt.plot(x,y)
plt.show()

plt.plot(x,amp)
plt.show()


