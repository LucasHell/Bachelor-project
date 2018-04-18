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


with open('transAmpl.txt', 'r') as inputFile:
	transitAmplitude = inputFile.readlines()

transitAmplitude = map(float, transitAmplitude)
for m in range(len(transitAmplitude)):
	if transitAmplitude[m] > 1 and transitAmplitude[m] < 66:
		amplCorr.append(transitAmplitude[m])
		
plt.hist(amplCorr,bins=11, rwidth=0.5)
plt.title("Histogram of the amplitudes of simulated TESS data")
plt.xlabel('Amplitude [min]')
plt.ylabel('#')
#~ plt.show()
plt.savefig('./plots/histo/ampl.png')
plt.clf()
		

with open('wtm-TESSTime.csv', 'r') as inputFile:
	data = inputFile.readlines()[1:]
	#~ RATess = [k.split(',')[0] for k in data]
	#~ decTess = [k.split(',')[1] for k in data]
	maxTime = [k.split(',')[2] for k in data]
	minTime = [k.split(',')[3] for k in data]
	medTime = [k.split(',')[4] for k in data]
	avgTime = [k.split(',')[5] for k in data]
	
with open('TESSTime.csv','r') as inputFile: 
	data = inputFile.readlines()
	RATess = [k.split(',')[0] for k in data]
	decTess = [k.split(',')[1] for k in data]
	
RATess = map(float, RATess)
decTess = map(float, decTess)
maxTime = map(float, maxTime)
minTime = map(float, minTime)
medTime = map(float, medTime)
transitAmplitude = map(float, transitAmplitude)


colorRange = np.linspace(float(np.min(transitAmplitude)), 100, float(np.max(transitAmplitude)))
#~ colorRange = np.linspace(0, 13, 13)


#~ norm = mp.colors.Normalize(
    #~ vmin=np.min(maxTime),
    #~ vmax=np.max(maxTime))
    
norm = mp.colors.Normalize(
    vmin=np.min(transitAmplitude),
    vmax=np.max(transitAmplitude))   
    
c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])



for i in range(len(RATess)):
	if decTess[i] > 0:
		RATess[i] = RATess[i] + 23
		decTess[i] = decTess[i] + 23
		if decTess[i] > 90:
			decTess[i] = -90 + (decTess[i]-90)
	else:
		RATess[i] = RATess[i] - 23
		decTess[i] = decTess[i] - 23
		if decTess[i] < -90:
			decTess[i] = 90 - (abs(decTess[i])-90)


for k in range(len(RATess)):
	if 300 <= RATess[k] <= 360 and 0 <= decTess[k] <= 30:
		print RATess[k]- 23, decTess[k] - 23, transitAmplitude[k] 
	elif 55 <= RATess[k] <= 70 and -35 <= decTess[k] <= -25:
		print RATess[k] + 23, decTess[k]+ 23, transitAmplitude[k] 

c = SkyCoord(ra=RATess*u.degree, dec=decTess*u.degree, frame='icrs')
ra_rad = c.ra.wrap_at(180 * u.deg).radian
dec_rad = c.dec.radian


plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m)
plt.scatter(ra_rad, dec_rad, s=7, c = transitAmplitude, cmap = cm.jet )
plt.savefig('plots/skymap_TESS_wrap')
plt.clf()

for i in range(len(transitAmplitude)):
	if float(transitAmplitude[i]) < 100:
		ra_cut.append(ra_rad[i])
		dec_cut.append(dec_rad[i])
		amp_cut.append(transitAmplitude[i])



norm = mp.colors.Normalize(
    vmin=np.min(amp_cut),
    vmax=np.max(amp_cut))   
    

c_m2 = mp.cm.cool
s_m2 = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m2.set_array([])

plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m2)
plt.scatter(ra_cut, dec_cut, s=7, c = amp_cut, cmap = cm.jet )
plt.savefig('plots/skymap_TESS_wrap_cutoff')
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
plt.savefig('./plots/histo/obsNumberMin.png')
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
	if float(transitAmplitude[i]) > 1 and float(transitAmplitude[i]) < 50:
		errorPlot.append(float(errorTiming[i])) 
		ampPlot.append(float(transitAmplitude[i])) 

	

x = np.linspace(1, np.amax(ampPlot))
y = x

plt.figure(figsize=(13, 13))
plt.scatter(np.log10(ampPlot), np.log10(errorPlot), s=20, c='black')
plt.xlabel('Amplitude [log(min)]', fontsize = 16)
plt.ylabel('Timing precision [log(min)]', fontsize = 16)
plt.title('Amplitude vs timing precision', fontsize = 18)
plt.plot(np.log10(x),np.log10(y))
plt.savefig('plots/ampErrorLog.png')
plt.clf()


plt.scatter(ampPlot, errorPlot, marker='d')
plt.xlabel('Amplitude [min]')
plt.ylabel('Timing precision [min]')
plt.title('Amplitude vs timing precision')
plt.savefig('plots/ampError.png')
plt.clf()

data = pd.read_table('ofir_table.txt', sep=';', skiprows=34, names=('KOI_num', 'newDetFlag', 'TTVfre', 'TTV+uncer', 'TTV-uncer', 'TTV_per', 'Delta_chi', 'chi_area', 'chi_single', 'chi_RMS', 'cho_correl', 'TTV_amp', 'TTV_amp+_uncer', 'TTV_amp-_uncer', 'TTV_ref', 'TTV_ref+_uncer', 'TTV_ref-_uncer', 'cofid', 'STD_error', '20', '21', '22', '23', '24'))
amp_ofir = data['TTV_amp']

#~ print amp_ofir
for i in range(len(amp_ofir)):
	if amp_ofir[i] > 5 and amp_ofir[i] < 66:
		ampOfirCorr.append(amp_ofir[i])

plt.hist(ampOfirCorr,bins=11, rwidth=0.5)
plt.title("Histogram of TTV amplitude of objects\nfrom the Ofir catalogue")
plt.xlabel('Amplitude [min]')
plt.ylabel('#')
plt.savefig('./plots/histo/ofir_amp.png')
plt.clf()	





	

