import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt

transitAmplitude = []
errorTiming = []
amp = []
err = []


with open('transAmplCheops.csv', 'r') as inputFile:
	data = inputFile.readlines()
	transitAmplitude = [k.split(',')[0] for k in data]
	
with open('ampErrorCheops.csv', 'r') as inputFile:
	data = inputFile.readlines()
	errorTiming = [k.split(',')[0] for k in data]


transitAmplitude = map(float, transitAmplitude)
errorTiming = map(float, errorTiming)



for l in range(len(transitAmplitude)):
	if transitAmplitude[l] > 0.0001:
		amp.append(transitAmplitude[l])
		err.append(errorTiming[l])
		
x = np.linspace(np.amin(amp), np.amax(amp))
y = x	
plt.scatter(np.log10(amp), np.log10(err))
plt.plot(np.log10(x),np.log10(y), c='black')
plt.savefig('plots/ampErrCheops.pdf')

