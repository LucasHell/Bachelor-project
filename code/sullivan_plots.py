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


with open('nasaDATA.csv','r') as inputFile: # read in data from csv file to respective arrays
	data = inputFile.readlines()[46:]
	kepID = [k.split(',')[1] for k in data]
	periodKep = [k.split(',')[5] for k in data]
	rPlanetKep = [k.split(',')[20] for k in data]
	name = [k.split(',')[2] for k in data]
	mStar = [k.split(',')[32] for k in data]
	numEpoch = [k.split(',')[8] for k in data]
	kepMag = [k.split(',')[37] for k in data]
	transitDur = [k.split(',')[14] for k in data]
	rStar = [k.split(',')[29] for k in data]
	
	
pTESS = 0
rTESS = 0
count1 = 0
count2 = 0

outputFile = open('TESSData.csv', 'w')
for i in range(len(periodSul)):
	for k in range(len(periodKep)):
		if rPlanetKep[k]:
			if periodSul[i] - periodRange <= float(periodKep[k]) <= periodSul[i] + periodRange and rPlanetSul[i] - radiusRange <= float(rPlanetKep[k]) <= rPlanetSul[i] + radiusRange:
				ratioP = float(periodKep[k]) / periodSul[i]
				ratioR = float(rPlanetKep[k]) / rPlanetSul[i]
				#~ j = i
				j = k
				count1 = 0
				count2 = 0
				
				while kepID[j] == kepID[j+1]:
					pTESS = float(periodKep[k+count1])*ratioP
					rTESS = float(rPlanetKep[k+count1])*ratioR
					if name[k] != '01' and count1 == 0:
						pNum = int(name[k][-2:])
						for l in range(1,pNum-1):
							periodTESS.append(float(periodKep[k-l])*ratioP)
							rPlanetTESS.append(float(rPlanetKep[k-l])*ratioR)
							outputFile.write(name[k-l] + ' ' +  kepID[k-l] + ' ' + str(pTESS) + ' ' + str(rTESS) + ' ' + mStar[k-l] + ' ' + numEpoch[k-l] + ' '  + transitDur[k-l] + ' ' + rStar[k-l] + ' ' + kepMag[k-l])
					
					periodTESS.append(pTESS)
					rPlanetTESS.append(rTESS)
					outputFile.write(name[k+count1] + ' ' +  kepID[k+count1] + ' ' + str(pTESS) + ' ' + str(rTESS) + ' ' + mStar[k+count1] + ' ' + numEpoch[k+count1] + ' '  + transitDur[k+count1] + ' ' + rStar[k+count1] + ' ' + kepMag[k+count1])
					print name[k+count1], kepID[k+count1], pTESS, rTESS, mStar[k+count1], numEpoch[k+count1], kepMag[k+count1], transitDur[k+count1], rStar[k+count1]
					count1 += 1
					j += 1			
				break
outputFile.close()
print len(rPlanetTESS), len(periodSul)					
periodTESS = np.log(periodTESS)
rPlanetTESS = np.log(rPlanetTESS)	
					
plt.scatter(rPlanetTESS,periodTESS,s=2)
plt.xlabel('Planet radius [$\log(R_{\oplus})$]', fontsize=12)
plt.ylabel('Period [$\log$ days]', fontsize=12)
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
	
