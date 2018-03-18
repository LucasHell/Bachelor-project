import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib  import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5


data = []
RATess = []			# Right ascension of TESS systems
decTess = []		# Declination of TESS systems
maxTime = []		# Max time that TESS observes the object
minTime = []		# Min time that TESS observes the object
medTime = []		# Med time that TESS observes the object
avgTime = []		# Average time that TESS observes the object
ampl = []			# amplitude
amplCorr = []		# amplitude where >0.1min is sorted out
errorPlot = []		
ampPlot = []
ra_rad = []
dec_rad = []


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

#~ plt.scatter( RATess, decTess, c = minTime, cmap = cm.jet )

norm = mp.colors.Normalize(
    vmin=np.min(maxTime),
    vmax=np.max(maxTime))
    
c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])



#~ c = SkyCoord(ra=RATess*u.degree, dec=decTess*u.degree, frame='icrs')
for i in range(len(RATess)):	
	c = SkyCoord(ra=RATess[i]*u.degree, dec=decTess[i]*u.degree, frame='icrs')
	c = c.galactic
	ra_rad.append(c.l)
	dec_rad.append(c.b)
	



#~ ra_rad = c.ra.wrap_at(180 * u.deg).radian
#~ dec_rad = c.dec.radian
#~ ra_rad = c_ecl.lon
#~ dec_rad = c_ecl.lat

plt.figure(figsize=(8,4.2))
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m)
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m)
plt.scatter(ra_rad, dec_rad, s=7,  c = maxTime, cmap = cm.jet )
plt.savefig('plots/skymap_TESS_wrap')
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
	
for i in range(len(errorTiming)):
	if transitAmplitude[i] == 'nan\n': 
		continue
	if float(transitAmplitude[i][:-2]) > 1 and float(transitAmplitude[i][:-2]) < 50:
		errorPlot.append(float(errorTiming[i][:-2])) 
		ampPlot.append(float(transitAmplitude[i][:-2])) 

	

x = np.linspace(1, np.amax(ampPlot))
y = x


plt.scatter(np.log10(ampPlot), np.log10(errorPlot), c='black')
plt.xlabel('Amplitude [log(min)]')
plt.ylabel('Timing precision [log(min)]')
plt.xlimit(2)
plt.ylimit(2)
plt.title('Amplitude vs timing precision')
plt.plot(np.log10(x),np.log10(y))
plt.savefig('plots/ampErrorLog.png')
plt.clf()

plt.scatter(ampPlot, errorPlot)
plt.xlabel('Amplitude [min]')
plt.ylabel('Timing precision [min]')
plt.title('Amplitude vs timing precision')
plt.savefig('plots/ampError.png')
plt.clf()

data = pd.read_table('sullivan_table.txt', sep=';', names=('KOI_num', 'newDetFlag', 'TTVfre', 'TTV+uncer', 'TTV-uncer', 'TTV_per', 'Delta_chi', 'chi_area', 'chi_single', 'chi_RMS', 'cho_correl', 'TTV_amp', 'TTV_amp+_uncer', 'TTV_amp-_uncer', 'TTV_ref', 'TTV_ref+_uncer', 'TTV_ref-_uncer', 'cofid', 'STD_error', '20', '21', '22', '23', '24'))
#~ print data
amp_ofir = data['KOI_num']
#~ print amp_ofir

#~ plt.hist(amp_ofir,bins=50)
#~ plt.title("Histogram of TTV amplitude of objects\nfrom the Ofir catalogue")
#~ plt.xlabel('Amplitude [min]')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/ofir_amp.png')
#~ plt.clf()	





	

