import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib  import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable


data = []
RATess = []			# Right ascension of TESS systems
decTess = []		# Declination of TESS systems
maxTime = []		# Max time that TESS observes the object
minTime = []		# Min time that TESS observes the object
medTime = []		# Med time that TESS observes the object
avgTime = []		# Average time that TESS observes the object
ampl = []			# amplitude
amplCorr = []		# amplitude where >0.1min is sorted out


with open('transAmpl.txt', 'r') as inputFile:
	ampl = inputFile.readlines()

ampl = map(float, ampl)
for m in range(len(ampl)):
	if ampl[m] > 0.5:
		amplCorr.append(ampl[m])
		
plt.hist(amplCorr,bins=10)
plt.title("Histogram of the amplitudes of simulated TESS data")
plt.xlabel('Amplitude [min]')
plt.ylabel('#')
#~ plt.show()
plt.savefig('./plots/histo/ampl.png')
plt.clf()
		

with open('wtm-TESSTime.csv', 'r') as inputFile:
	data = inputFile.readlines()[1:]
	RATess = [k.split(',')[0] for k in data]
	decTess = [k.split(',')[1] for k in data]
	maxTime = [k.split(',')[2] for k in data]
	minTime = [k.split(',')[3] for k in data]
	medTime = [k.split(',')[4] for k in data]
	avgTime = [k.split(',')[5] for k in data]
	
RATess = map(float, RATess)
decTess = map(float, decTess)
maxTime = map(float, maxTime)
minTime = map(float, minTime)
medTime = map(float, medTime)

colorRange = np.linspace(0, 13, 13)

plt.scatter( RATess, decTess, c = minTime, cmap = cm.jet )

norm = mp.colors.Normalize(
    vmin=np.min(minTime),
    vmax=np.max(minTime))
    
c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])

plt.colorbar(s_m)
plt.xlabel('Right Ascension [$^{\circ}$]', fontsize=12)
plt.ylabel('Declination [$^{\circ}$]', fontsize=12)
plt.savefig('./plots/RA_Dec_Min.png')
plt.clf()

plt.hist(minTime,bins=50)
plt.title("Histogram of times that each object are observed by TESS")
plt.xlabel('# of times observed')
plt.ylabel('#')
plt.savefig('./plots/histo/obsNumberMin.png')
plt.clf()

	
data = pd.read_table('sullivan_table.txt', delim_whitespace=True, names=('RA', 'Dec', 'Rp', 'P', 'P insu', 'Rad v', 'Rs', 'Teff', 'Vmag', 'ICmag', 'Jmag', 'KSmag', 'Dmod', 'Dil p', 'devi flux', 'Sig-noi', 'NumPl'))

effTemp = data['Teff']
sRad = data['Rs']
vMag = data['Vmag']
ICMag = data['ICmag']
jMag = data['Jmag']
KSMag = data['KSmag']


			

plt.hist(effTemp,bins=50)
plt.title("Histogram of the effective temperatures of objects\nfrom the Sullivan catalogue")
plt.xlabel('Effective temperature [K$^\circ$]')
plt.ylabel('#')
plt.savefig('./plots/histo/effTemp.png')
plt.clf()

plt.hist(sRad,bins=50)
plt.title("Histogram of the star radii of objects\nfrom the Sullivan catalogue")
plt.xlabel('Star radius [R$_{\odot}$]')
plt.ylabel('#')
plt.savefig('./plots/histo/sRad.png')
plt.clf()

plt.hist(vMag,bins=50)
plt.title("Histogram of the V-band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('V-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/vMag.png')
plt.clf()

plt.hist(ICMag,bins=50)
plt.title("Histogram of the I_C_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('IC-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/ICMag.png')
plt.clf()

plt.hist(jMag,bins=50)
plt.title("Histogram of the J-band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('J-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/jMag.png')
plt.clf()

plt.hist(KSMag,bins=50)
plt.title("Histogram of the K_s_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('KS-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/KSMag.png')
plt.clf()

with open('ampError.txt','r') as inputFile:
	errorTiming = inputFile.readlines()
	
with open('transAmpl.txt','r') as inputFile:
	transitAmplitude = inputFile.readlines()

print errorTiming[0], transitAmplitude[0]
plt.scatter(float(transitAmplitude), float(errorTiming))
plt.xlabel('Amplitude [min]')
plt.ylabel('Timing precision [min]')
plt.title('Amplitude vs timing precision')
plt.savefig('plots/ampError.png')



	

