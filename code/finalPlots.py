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
from numpy import median


data = []
RATess = []					# Right ascension of TESS systems
decTess = []				# Declination of TESS systems
maxTime = []				# Max time that TESS observes the object
minTime = []				# Min time that TESS observes the object
medTime = []				# Med time that TESS observes the object
avgTime = []				# Average time that TESS observes the object
transitAmplitude = []		# amplitude
amplCorr = []				# amplitude where <0.1min is sorted out
errorPlot = []		
ampPlot = []
ra_rad = []
dec_rad = []
ra_radTime = []
dec_radTime = []
ampOfirCorr = []
ra_cut = []
dec_cut = []
amp_cut = []
timeCut = []
beta = []
lamb = []
numP = []
sysName = []
ampl = []
periodFrac = []
amplFil = []
periodFracFil = []
amplDouble = []
periodFracDouble = []
amplFilDouble = []
periodFracFilDouble = []
rPlanet = []
betaTime = []
lambTime = []



with open('transAmpl.csv', 'r') as inputFile:
	data = inputFile.readlines()
	transitAmplitude = [k.split(',')[0] for k in data]
	#~ RATess = [k.split(',')[1] for k in data]
	#~ decTess = [k.split(',')[2] for k in data]
	rPlanet = [k.split(',')[4] for k in data]
	#~ sysName = [k.split(',')[3] for k in data]


transitAmplitude = map(float, transitAmplitude)
for m in range(len(transitAmplitude)):
	if transitAmplitude[m] > 5 and transitAmplitude[m] < 300:
		amplCorr.append(transitAmplitude[m])
		
		

y,binEdges = np.histogram(amplCorr,bins=11)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
menStd     = np.sqrt(y)
width = 5
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

with open('wtm-RA_dec_sys.csv', 'r') as inputFile:
	data = inputFile.readlines()[1:]
	#~ RATessTime = [k.split(',')[0] for k in data]
	#~ decTessTime = [k.split(',')[1] for k in data]
	maxTime = [k.split(',')[2] for k in data]
	minTime = [k.split(',')[3] for k in data]
	medTime = [k.split(',')[4] for k in data]
	avgTime = [k.split(',')[5] for k in data]
	
with open('RA_dec_p.csv','r') as inputFile: 
	data = inputFile.readlines()
	RATess = [k.split(',')[0] for k in data]
	decTess = [k.split(',')[1] for k in data]

with open('RA_dec_sys.csv','r') as inputFile: 
	data = inputFile.readlines()
	RATessTime = [k.split(',')[0] for k in data]
	decTessTime = [k.split(',')[1] for k in data]
	
RATess = map(float, RATess)
decTess = map(float, decTess)
RATessTime = map(float, RATessTime)
decTessTime = map(float, decTessTime)
maxTime = map(float, maxTime)
minTime = map(float, minTime)
medTime = map(float, medTime)
transitAmplitude = map(float, transitAmplitude)





coordsEcl = open('coordsEcl.csv', 'w')
coordsEcl.write('RA' + ',' + 'dec' + ',' + 'long' + ',' + 'lat' + ',' + 'Name' + '\n')
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
	coordsEcl.write(repr(math.degrees(RATess[i])) + ',' + repr(math.degrees(decTess[i])) + ',' + repr(lamb[i]) + ',' + repr(beta[i]) + '\n')  

c = SkyCoord(lon=lamb*u.deg, lat=beta*u.deg, frame='heliocentrictrueecliptic')
ra_rad = c.lon.wrap_at(180 * u.deg).rad			
dec_rad = c.lat.rad


for i in range(len(RATessTime)):
	RATessTime[i] = math.radians(RATessTime[i])
	decTessTime[i] = math.radians(decTessTime[i])


	betaTime.append(math.asin(math.cos(tilt)*math.sin(decTessTime[i]) - math.sin(RATessTime[i])*math.cos(decTessTime[i])*math.sin(tilt)))
	lamb1 = (math.sin(tilt)*math.sin(decTessTime[i]) + math.sin(RATessTime[i])*math.cos(decTessTime[i])*math.cos(tilt))/math.cos(beta[i])
	lamb2 = (math.cos(RATessTime[i])*math.cos(decTessTime[i]))/math.cos(beta[i])
	lambTime.append(math.degrees(math.atan2(lamb2, lamb1)))
	betaTime[i] = math.degrees(betaTime[i])
	coordsEcl.write(repr(math.degrees(RATessTime[i])) + ',' + repr(math.degrees(decTessTime[i])) + ',' + repr(lambTime[i]) + ',' + repr(betaTime[i]) + '\n')  

c = SkyCoord(lon=lambTime*u.deg, lat=betaTime*u.deg, frame='heliocentrictrueecliptic')
ra_radTime = c.lon.wrap_at(180 * u.deg).rad			
dec_radTime = c.lat.rad

norm = mp.colors.Normalize(
    vmin=np.min(maxTime),
    vmax=np.max(maxTime))
    
#~ norm = mp.colors.Normalize(
    #~ vmin=np.min(transitAmplitude),
    #~ vmax=np.max(transitAmplitude)) 

c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.jet, norm=norm)
s_m.set_array([])

#~ print len(ra_radTime), len(dec_radTime), len(maxTime)
plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.colorbar(s_m)
plt.axhspan(math.radians(-40), math.radians(40), facecolor='g', alpha=0.1)
plt.scatter(ra_radTime, dec_radTime, s=7, c = maxTime, cmap = cm.jet )
textstr = 'Number of observations'
plt.figtext(0.88, 0.7, textstr, fontsize=12, rotation=90)
plt.savefig('plots/skymap_TESS_numObs.pdf')
plt.clf()
countAmp = 0
#~ print transitAmplitude.index(np.amax(transitAmplitude))


#~ print len(ra_rad), len(dec_rad), len(transitAmplitude), len(maxTime)
#~ print ra_rad[1]
for i in range(len(transitAmplitude)):
	if float(transitAmplitude[i]) < 3000:
		ra_cut.append(ra_rad[i])
		dec_cut.append(dec_rad[i])
		amp_cut.append(transitAmplitude[i])
		#~ if -90 < math.degrees(ra_rad[i]) < -60 and 15 < math.degrees(dec_rad[i]) < 45:
			#~ print i, transitAmplitude[i]
		#~ timeCut.append(maxTime[i])
		#~ if -40 < dec_rad[i] < 40 and float(transitAmplitude[i]) > 20 and dec_rad[i] != dec_rad[i-1]:
			#~ countAmp += 1
			#~ print math.degrees(ra_rad[i]), math.degrees(dec_rad[i]), np.log10(float(transitAmplitude[i]))


#~ print countAmp

for l in range(len(amp_cut)):
	if amp_cut[l] < 0.0001:
		amp_cut[l] = -5
	else:
		amp_cut[l] = np.log10(amp_cut[l])
    
#~ norm = mp.colors.Normalize(
    #~ vmin=np.min(timeCut),
    #~ vmax=np.max(timeCut))   


norm = mp.colors.Normalize(
    vmin=np.min(-1),
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
textstr = 'TTV amplitude [min]'
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
	
	


countSys = 0
#~ print len(errorTiming), len(transitAmplitude)
for i in range(len(errorTiming)):
	if transitAmplitude[i] == 'nan\n' or transitAmplitude[i] == '0\n': 
		#~ countSys += 1
		continue
	elif float(transitAmplitude[i]) < 0.001:
		errorPlot.append(float(errorTiming[i])) 
		ampPlot.append(np.random.normal(0.2, 0.05))
	else:
		errorPlot.append(float(errorTiming[i])) 
		ampPlot.append(float(transitAmplitude[i])) 
		
		
		
#~ print np.log10(median(errorPlot))
#~ print countSys
x = np.linspace(np.amin(ampPlot), np.amax(ampPlot))
y = x

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(13, 13))
plt.scatter(np.log10(ampPlot), np.log10(errorPlot), s=20)
plt.xlabel('log$_{10}$[Amplitude (min)]', fontsize = 20)
plt.ylabel('log$_{10}$[Error (min)]', fontsize = 20)
#~ plt.title('Amplitude vs Error', fontsize = 18)
plt.plot(np.log10(x),np.log10(y), c='black')
plt.savefig('plots/ampErrorLog.pdf')
plt.clf()


countAmp = 0
countNo = 0
for o in range(len(ampPlot)):
	if ampPlot[o] > errorPlot[o] and 75 < math.degrees(dec_cut[o]) < 90 or -75 < math.degrees(dec_cut[o]) < -90:
		countAmp += 1
		countNo += 1
	elif 75 < math.degrees(dec_cut[o]) < 90 or -75 < math.degrees(dec_cut[o]) < -90:
		countNo += 1
print math.degrees(dec_cut[0]), math.degrees(dec_cut[1])
print "Count: ", countAmp, countNo, len(dec_cut)

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

#~ with open('AmplPeriod.csv','r') as inputFile: 
	#~ data = inputFile.readlines()
	#~ periodFrac = [k.split(',')[0] for k in data]
	#~ ampl = [k.split(',')[1] for k in data]
	#~ numP = [k.split(',')[2] for k in data]
	
#~ periodFrac = map(float, periodFrac)
#~ ampl = map(float, ampl)
#~ numP = map(float, numP)

#~ for i in range(len(ampl)):
	#~ if ampl[i] > 0.0001:
		#~ amplFil.append(ampl[i])
		#~ periodFracFil.append(periodFrac[i])
		
#~ plt.scatter(periodFracFil, np.log10(amplFil))
#~ plt.xlabel('Period ratio')
#~ plt.ylabel('Amplitude [min]')
#~ plt.savefig('./plots/ampPrat.pdf')
#~ plt.clf()

#~ plt.scatter(periodFrac, numP)
#~ plt.xlabel('Period ratio')
#~ plt.ylabel('Multiplicity')
#~ plt.savefig('./plots/multiPrat.pdf')
#~ plt.clf()



with open('AmplPeriodDouble.csv','r') as inputFile: 
	data = inputFile.readlines()
	periodFracDouble = [k.split(',')[0] for k in data]
	amplDouble = [k.split(',')[1] for k in data]
	rPlanetDouble = [k.split(',')[2] for k in data]



periodFracDouble = map(float, periodFracDouble)
amplDouble = map(float, amplDouble)
rPlanetDouble = map(float, rPlanetDouble)

plt.hist(periodFracDouble, bins=20)
plt.xlabel('Period ratio')
plt.ylabel('Number of planets')
plt.savefig('./plots/PratHisto.pdf')
plt.clf()

#~ plt.scatter(np.log10(rPlanetDouble), np.log10(amplDouble))
#~ plt.xlabel('log$_{10}$[Radius (R$_{\oplus}$)]')
#~ plt.ylabel('log$_{10}$[Amplitude (min)]')
#~ plt.savefig('./plots/ampRadiusDouble.pdf')
#~ plt.clf()

#~ for i in range(len(amplDouble)):
	#~ if amplDouble[i] > 0.0001 and periodFracDouble[i] < 3:
		#~ amplFilDouble.append(amplDouble[i])
		#~ periodFracFilDouble.append(periodFracDouble[i])
		
#~ plt.scatter(periodFracFilDouble, np.log10(amplFilDouble))
#~ plt.xlabel('Period ratio')
#~ plt.ylabel('Amplitude [min]')
#~ plt.savefig('./plots/ampPratDouble.pdf')
#~ plt.clf()

rPlanet = map(float, rPlanet)
transitAmplitude = map(float, transitAmplitude)
rPlanet2 = []
transAmp2 = []
for i in range(len(rPlanet)):
	if 0.0001 < transitAmplitude[i] < 1000:
		rPlanet2.append(rPlanet[i])
		transAmp2.append(transitAmplitude[i])
#~ transitAmplitude = [i for i in transitAmplitude if i < 1000]

#~ print np.amax(rPlanet2), np.amax(transAmp2)
plt.scatter(np.log10(rPlanet2), np.log10(transAmp2))
plt.xlabel('log$_{10}$[Radius (R$_{\oplus}$)]')
plt.ylabel('log$_{10}$[Amplitude (min)]')
plt.savefig('./plots/ampRadius.pdf')
plt.clf()



	

