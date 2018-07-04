import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib  import cm
from numpy import unravel_index

transitAmplitude = []
errorTiming = []
amp = []
err = []
periodFrac = []
rPlanet = []
ampl = []
periodFracCut = []
rPlanetCut = []
amplCut = []
periodFracHist = []
rPlanetHist = []
rPlanetHist2 = []
periodFracHist2 = []
H = []




with open('transAmplCheops.csv', 'r') as inputFile:
	data = inputFile.readlines()
	transitAmplitude = [k.split(',')[0] for k in data]

with open('AmplPeriod.csv', 'r') as inputFile:
	data = inputFile.readlines()
	periodFrac = [k.split(',')[0] for k in data]
	rPlanet = [k.split(',')[1] for k in data]
	ampl = [k.split(',')[2] for k in data]
	error = [k.split(',')[3] for k in data]
	
	
with open('ampErrorCheops.csv', 'r') as inputFile:
	data = inputFile.readlines()
	errorTiming = [k.split(',')[0] for k in data]


transitAmplitude = map(float, transitAmplitude)
errorTiming = map(float, errorTiming)
periodFrac = map(float, periodFrac)
rPlanet = map(float, rPlanet)
ampl = map(float, ampl)


for l in range(len(errorTiming)):
	if transitAmplitude[l] > 0.0001:
		amp.append(transitAmplitude[l])
		err.append(errorTiming[l]*60)
		
x = np.linspace(0.1, 100)
y = x	
plt.scatter(np.log10(amp), np.log10(err))
plt.plot(np.log10(x),np.log10(y), c='black')
plt.xlabel('log$_{10}$[Amplitude (min)]')
plt.ylabel('log$_{10}$[Error (min)]')
plt.savefig('plots/ampErrCheops.pdf')
plt.clf()


for l in range(len(rPlanet)):
	if rPlanet[l] < 20 and 0.003 < ampl[l] < 650:
		periodFracCut.append(periodFrac[l])
		rPlanetCut.append(rPlanet[l])
		if ampl[l] < 10^(-10) or str(ampl[l]) == 'nan' or ampl[l] < 0:
			amplCut.append(10^(-10))
		else:
			amplCut.append(np.log10(ampl[l]))
		#~ print np.amin(amplCut), ampl[l]
#~ print amplCut
#~ print np.amin(ampl)

norm = mp.colors.Normalize(
    vmin=np.amin(amplCut),
    vmax=np.amax(amplCut))


c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])


plt.scatter(periodFracCut, rPlanetCut, c = amplCut, cmap = cm.jet, alpha = 0.6)
plt.colorbar(s_m)
plt.xlabel('Period ratio')
plt.ylabel('Planet radius [R$_{\oplus}$]')
textstr = 'log$_{10}$[TTV amplitude (min)]'
plt.figtext(0.88, 0.7, textstr, fontsize=12, rotation=90)
plt.savefig('plots/perFracRad.pdf')
plt.clf()

norm = mp.colors.Normalize(
    vmin=np.amin(amplCut),
    vmax=np.amax(amplCut))
c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])



for i in range(len(ampl)):
	if ampl[i] > float(error[i])*60 and rPlanet[i] < 20:
		periodFracHist.append(periodFrac[i])
		rPlanetHist.append(rPlanet[i])
	if rPlanet[i] < 20:
		periodFracHist2.append(periodFrac[i])
		rPlanetHist2.append(rPlanet[i])

print len(periodFracHist), len(rPlanetHist2)
HAll, xedges, yedges = np.histogram2d(periodFracHist2, rPlanetHist2, bins = 12)
HErr, _,_ = np.histogram2d(periodFracHist, rPlanetHist, bins = 12)
extent = [1,6, yedges[0], yedges[-1]]
#~ print HErr[2,1]
#~ print HAll[2,1]

for i in range(0, 12):
	H.append([])
	for l in range(0, 12):
		if HAll[i,l] == 0:
			HAll[i,l] = 1
		if HErr[i,l] > HAll[i,l]:
			H[i].append(0.5)
		else:
			H[i].append(HErr[i,l]/HAll[i,l]*100)

#~ print H[2][1]
#~ print HAll
#~ print 
#~ print HErr
plt.imshow(H, extent=extent, interpolation='nearest', aspect = 'auto')
plt.xlabel('Period ratio', fontsize = 14)
plt.ylabel('Planet radius [R$_{\oplus}$]', fontsize = 14)
plt.colorbar()
textstr = 'TTV probability [%]'
plt.figtext(0.88, 0.7, textstr, fontsize=12, rotation=90)
plt.savefig('plots/2dhist3.pdf')
plt.clf()


plt.hist2d(periodFracHist, rPlanetHist, bins = 10)
plt.colorbar()
plt.savefig('plots/2dhistErr.png')
plt.clf()


plt.hist2d(periodFracHist2, rPlanetHist2, bins = 10)
plt.colorbar()
plt.savefig('plots/2dhistAll.png')
plt.clf()


plt.scatter(periodFracCut, amplCut)
plt.gca().set_xlim(right = 3)
plt.savefig('plots/perAml.pdf')
plt.clf()


























