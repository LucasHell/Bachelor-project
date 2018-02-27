import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import astropy.table as Table

effTemp = []		# effective temperature
sRad = [] 			# stellar radius
vMag = []			# V-band magnitude (all magnitudes are apparent magnitude)
ICMag = []			# I_C magnitude
jMag = []			# J-band magnitude
KSMag = []			# K_S magnitude
periodSul = []		# period from Sullivan catalogue
rPlanetSul = []		# planet radius from Sullivan catalogue
periodTESS = []		# period calculated from Kepler and Sullivan data
rPlanetTESS = []	# planet radius calculated from Kepelr and Sullivan data
kepID = []			# Kepler ID of planets

data = pd.read_table('sullivan_table.txt', delim_whitespace=True, names=('RA', 'Dec', 'Rp', 'P', 'P insu', 'Rad v', 'Rs', 'Teff', 'Vmag', 'ICmag', 'Jmag', 'KSmag', 'Dmod', 'Dil p', 'devi flux', 'Sig-noi', 'NumPl'))

effTemp = data['Teff']
sRad = data['Rs']
vMag = data['Vmag']
ICMag = data['ICmag']
jMag = data['Jmag']
KSMag = data['KSmag']
periodSul = data['P']
rPlanetSul = data['Rp']

periodRange = 1
radiusRange = 1

count = 0

with open('nasaDATA.csv','r') as inputFile: # read in data from csv file to respective arrays
	data = inputFile.readlines()[46:]
	kepID = [k.split(',')[1] for k in data]
	periodKep = [k.split(',')[5] for k in data]
	rPlanetKep = [k.split(',')[20] for k in data]
	
	
	
#~ with open('numberPlanets.csv','r') as inputFile:
	#~ data = inputFile.readlines()[46:]
	#~ numP = [k.split(' ')[0] for k in data]



for i in range(len(periodSul)):
	for k in range(len(periodKep)):
		if rPlanetKep[k]:
			if periodSul[i] - periodRange <= float(periodKep[k]) <= periodSul[i] + periodRange and rPlanetSul[i] - radiusRange <= float(rPlanetKep[k]) <= rPlanetSul[i] + radiusRange:
				ratioP = float(periodKep[k]) / periodSul[i]
				ratioR = float(rPlanetKep[k]) / rPlanetSul[i]
				#~ print i, k
				#~ print i, j
				j = i
				count = 0
				while kepID[j] == kepID[j+1]:
					periodTESS.append(float(periodKep[k+count])*ratioP)
					rPlanetTESS.append(float(rPlanetKep[k+count])*ratioR)
					count += 1
					j += 1
					print i, k
				break
			
print len(rPlanetTESS)						
periodTESS = np.log(periodTESS)
rPlanetTESS = np.log(rPlanetTESS)	
print len(rPlanetTESS)		
					
plt.scatter(rPlanetTESS,periodTESS)
plt.show()			

#~ plt.hist(effTemp,bins=50)
#~ plt.title("Histogram of the effective temperatures of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('Effective temperature [K$^\circ$]')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/effTemp.png')
#~ plt.clf()

#~ plt.hist(sRad,bins=50)
#~ plt.title("Histogram of the star radii of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('Star radius [R$_{\odot}$]')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/sRad.png')
#~ plt.clf()

#~ plt.hist(vMag,bins=50)
#~ plt.title("Histogram of the V-band magnitude of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('V-band apparent magnitude')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/vMag.png')
#~ plt.clf()

#~ plt.hist(ICMag,bins=50)
#~ plt.title("Histogram of the I_C_ band magnitude of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('IC-band apparent magnitude')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/ICMag.png')
#~ plt.clf()

#~ plt.hist(jMag,bins=50)
#~ plt.title("Histogram of the J-band magnitude of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('J-band apparent magnitude')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/jMag.png')
#~ plt.clf()

#~ plt.hist(KSMag,bins=50)
#~ plt.title("Histogram of the K_s_ band magnitude of objects\nfrom the Sullivan catalogue")
#~ plt.xlabel('KS-band apparent magnitude')
#~ plt.ylabel('#')
#~ plt.savefig('./plots/histo/KSMag.png')
	
