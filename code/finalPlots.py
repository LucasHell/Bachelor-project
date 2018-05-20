import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib  import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import matplotlib.lines as mlines
import math


data = []
RATess = []					# Right ascension of TESS systems
decTess = []				# Declination of TESS systems
maxTime = []				# Max time that TESS observes the object
minTime = []				# Min time that TESS observes the object
medTime = []				# Med time that TESS observes the object
avgTime = []				# Average time that TESS observes the object
transitAmplitude = []		# amplitude
amplCorr = []				# amplitude where >0.1min is sorted out
errorPlot = []		
ampPlot = []
ra_rad = []
dec_rad = []
ampOfirCorr = []
ra_cut = []
dec_cut = []
amp_cut = []
timeCut = []
beta = []
lamb = []



with open('transAmpl.csv', 'r') as inputFile:
	data = inputFile.readlines()
	transitAmplitude = [k.split(',')[0] for k in data]
	RATess = [k.split(',')[1] for k in data]
	decTess = [k.split(',')[2] for k in data]

transitAmplitude = map(float, transitAmplitude)
for m in range(len(transitAmplitude)):
	if transitAmplitude[m] > 5 and transitAmplitude[m] < 400:
		amplCorr.append(transitAmplitude[m])
		
		

y,binEdges = np.histogram(amplCorr,bins=20)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
menStd     = np.sqrt(y)
width = 10
plt.bar(bincenters, y, width=width, yerr=menStd, error_kw=dict(ecolor='black', lw=1, capsize=4, capthick=1))
#~ plt.title("Histogram of the amplitudes of simulated Ofir data")
plt.xlabel('Amplitude [min]')
plt.ylabel('Number of planets')
#~ plt.show()
plt.savefig('./plots/histo/ampl.pdf')
plt.clf()
		
transSort = sorted(transitAmplitude, key=float, reverse=False)
outTransFile = open('transAmpSort.csv', 'w')
for i in range(len(transSort)):
	outTransFile.write(str(transSort[i]) + '\n')
outTransFile.close()
with open('wtm-TESSTime.csv', 'r') as inputFile:
	data = inputFile.readlines()[1:]
	#~ RATess = [k.split(',')[0] for k in data]
	#~ decTess = [k.split(',')[1] for k in data]
	maxTime = [k.split(',')[2] for k in data]
	minTime = [k.split(',')[3] for k in data]
	medTime = [k.split(',')[4] for k in data]
	avgTime = [k.split(',')[5] for k in data]
	
#~ with open('RA_dec_p.csv','r') as inputFile: 
	#~ data = inputFile.readlines()
	#~ RATess = [k.split(',')[0] for k in data]
	#~ decTess = [k.split(',')[1] for k in data]
	
RATess = map(float, RATess)
decTess = map(float, decTess)
maxTime = map(float, maxTime)
minTime = map(float, minTime)
medTime = map(float, medTime)
transitAmplitude = map(float, transitAmplitude)






tilt = math.radians(23.439281)
lamb1 = 0
lamb2 = 0
for i in range(len(RATess)):
	RATess[i] = math.radians(RATess[i])
	decTess[i] = math.radians(decTess[i])


	beta.append(math.asin(math.cos(tilt)*math.sin(decTess[i]) - math.sin(RATess[i])*math.cos(decTess[i])*math.sin(tilt)))
	lamb1 = (math.sin(tilt)*math.sin(decTess[i]) + math.sin(RATess[i])*math.cos(decTess[i])*math.cos(tilt))/math.cos(beta[i])
	lamb2 = (math.cos(RATess[i])*math.cos(decTess[i]))/math.cos(beta[i])
	lamb.append(math.degrees(math.atan2(lamb2, lamb1)))
	beta[i] = math.degrees(beta[i])

c = SkyCoord(lon=lamb*u.deg, lat=beta*u.deg, frame='heliocentrictrueecliptic')
ra_rad = c.lon.wrap_at(180 * u.deg).rad			
dec_rad = c.lat.rad



norm = mp.colors.Normalize(
    vmin=np.min(maxTime),
    vmax=np.max(maxTime))
    
norm = mp.colors.Normalize(
    vmin=np.min(transitAmplitude),
    vmax=np.max(transitAmplitude)) 

c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])

#~ print len(ra_rad), len(dec_rad), len(transitAmplitude)
plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m)
#~ plt.scatter(ra_rad, dec_rad, s=7, c = transitAmplitude, cmap = cm.jet )
plt.savefig('plots/skymap_TESS_wrap')
plt.clf()
countAmp = 0
print transitAmplitude.index(np.amax(transitAmplitude))
#~ print transitAmplitude.index(np.amax(transitAmplitude))
#~ print math.degrees(RATess[462]), math.degrees(decTess[462])
print len(ra_rad), len(dec_rad), len(transitAmplitude)
for i in range(len(transitAmplitude)):
	if float(transitAmplitude[i]) < 3000:
		ra_cut.append(ra_rad[i] * -1)
		dec_cut.append(dec_rad[i])
		amp_cut.append(transitAmplitude[i])
		timeCut.append(maxTime[i])
		if -40 < dec_rad[i] < 40 and float(transitAmplitude[i]) > 20 and dec_rad[i] != dec_rad[i-1]:
			countAmp += 1
			print math.degrees(ra_rad[i]), math.degrees(dec_rad[i]), np.log10(float(transitAmplitude[i]))


print countAmp

for l in range(len(amp_cut)):
	if amp_cut[l] < 0.0001:
		amp_cut[l] = 0.0001
	else:
		amp_cut[l] = np.log10(amp_cut[l])
    
#~ norm = mp.colors.Normalize(
    #~ vmin=np.min(timeCut),
    #~ vmax=np.max(timeCut))   


norm = mp.colors.Normalize(
    vmin=np.min(amp_cut),
    vmax=np.max(amp_cut))  
    
   


c_m2 = mp.cm.cool
s_m2 = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m2.set_array([])

plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
#~ plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m2)
plt.axhspan(math.radians(-40), math.radians(40), facecolor='g', alpha=0.1)
plt.scatter(ra_cut, dec_cut, s=7, c = amp_cut, cmap = cm.jet, alpha = 0.5)
textstr = 'Number of observations'
plt.figtext(0.88, 0.7, textstr, fontsize=12, rotation=90)
plt.savefig('plots/skymap_TESS_amp.pdf')
plt.clf()



#~ plt.colorbar(s_m)
#~ plt.xlabel('Right Ascension [$^{\circ}$]', fontsize=12)
#~ plt.ylabel('Declination [$^{\circ}$]', fontsize=12)
#~ plt.savefig('./plots/RA_Dec_Min.png')
#~ plt.clf()

plt.hist(minTime,bins=50)
plt.title("Histogram of times that each object are observed by TESS")
plt.xlabel('# of times observed')
plt.ylabel('#')
plt.savefig('./plots/histo/obsNumberMin.pdf')
plt.clf()

	
data = pd.read_table('sullivan_table.txt', delim_whitespace=True, names=('RA', 'Dec', 'Rp', 'P', 'P insu', 'Rad v', 'Rs', 'Teff', 'Vmag', 'ICmag', 'Jmag', 'KSmag', 'Dmod', 'Dil p', 'devi flux', 'Sig-noi', 'NumPl'))

effTemp = data['Teff']
sRad = data['Rs']
vMag = data['Vmag']
ICMag = data['ICmag']
jMag = data['Jmag']
KSMag = data['KSmag']


			

plt.hist(effTemp,bins=50, rwidth=0.5)
plt.title("Histogram of the effective temperatures of objects\nfrom the Sullivan catalogue")
plt.xlabel('Effective temperature [K$^\circ$]')
plt.ylabel('#')
plt.savefig('./plots/histo/effTemp.png')
plt.clf()

plt.hist(sRad,bins=50, rwidth=0.5)
plt.title("Histogram of the star radii of objects\nfrom the Sullivan catalogue")
plt.xlabel('Star radius [R$_{\odot}$]')
plt.ylabel('#')
plt.savefig('./plots/histo/sRad.png')
plt.clf()

plt.hist(vMag,bins=50, rwidth=0.5)
plt.title("Histogram of the V-band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('V-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/vMag.png')
plt.clf()

plt.hist(ICMag,bins=50, rwidth=0.5)
plt.title("Histogram of the I_C_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('IC-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/ICMag.png')
plt.clf()

plt.hist(jMag,bins=50, rwidth=0.5)
plt.title("Histogram of the J-band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('J-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/jMag.png')
plt.clf()

plt.hist(KSMag,bins=50, rwidth=0.5)
plt.title("Histogram of the K_s_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.xlabel('KS-band apparent magnitude')
plt.ylabel('#')
plt.savefig('./plots/histo/KSMag.png')
plt.clf()

with open('ampError.txt','r') as inputFile:
	errorTiming = inputFile.readlines()


for i in range(len(errorTiming)):
	if transitAmplitude[i] == 'nan\n' or transitAmplitude[i] == '0\n': 
		continue
	if float(transitAmplitude[i]) > 0.001 and float(transitAmplitude[i]) < 200:
		errorPlot.append(float(errorTiming[i])) 
		ampPlot.append(float(transitAmplitude[i])) 

	

x = np.linspace(np.amin(ampPlot), np.amax(ampPlot))
y = x

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(13, 13))
plt.scatter(np.log10(ampPlot), np.log10(errorPlot), s=20)
plt.xlabel('log$_{10}$[Amplitude (min)]', fontsize = 18)
plt.ylabel('log$_{10}$[Error (min)]', fontsize = 18)
#~ plt.title('Amplitude vs Error', fontsize = 18)
plt.plot(np.log10(x),np.log10(y), c='black')
plt.savefig('plots/ampErrorLog.pdf')
plt.clf()


plt.scatter(ampPlot, errorPlot, marker='d')
plt.xlabel('Amplitude [min]')
plt.ylabel('Timing precision [min]')
#~ plt.title('Amplitude vs timing precision')
plt.savefig('plots/ampError.pdf')
plt.clf()

data = pd.read_table('ofir_table.txt', sep=';', skiprows=34, names=('KOI_num', 'newDetFlag', 'TTVfre', 'TTV+uncer', 'TTV-uncer', 'TTV_per', 'Delta_chi', 'chi_area', 'chi_single', 'chi_RMS', 'cho_correl', 'TTV_amp', 'TTV_amp+_uncer', 'TTV_amp-_uncer', 'TTV_ref', 'TTV_ref+_uncer', 'TTV_ref-_uncer', 'cofid', 'STD_error', '20', '21', '22', '23', '24'))
amp_ofir = data['TTV_amp']

for i in range(len(amp_ofir)):
	if amp_ofir[i] > 1 and amp_ofir[i] < 66:
		ampOfirCorr.append(amp_ofir[i])

plt.hist(ampOfirCorr,bins=11, rwidth=0.5)
#~ plt.title("Histogram of TTV amplitude of objects\nfrom the Ofir catalogue")
plt.xlabel('Amplitude [min]')
plt.ylabel('#')
plt.savefig('./plots/histo/ofir_amp.pdf')
plt.clf()	








	

