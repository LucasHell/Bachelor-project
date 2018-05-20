import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mp
from matplotlib  import cm
import matplotlib.lines as mlines

effTemp = []		# effective temperature
sRad = [] 			# stellar radius
vMag = []			# V-band magnitude (all magnitudes are apparent magnitude)
ICMag = []			# I_C magnitude
jMag = []			# J-band magnitude
KSMag = []			# K_S magnitude
periodSul = []		# period from Sullivan catalogue
rPlanetSul = []		# planet radius from Sullivan catalogue
periodTESS = []		# period calculated from Kepler and Sullivan data
rPlanetTESS = []	# planet radius calculated from Kepler and Sullivan data
kepID = []			# Kepler ID of planets
RA = []				# Right ascension
dec = []			# Declination
eRay = []			# Eccentricity (Rayleigh distribution)
omegaUni = []		# Omega (uniform distribution)
numP = []			# Number of planets
pNum = 0			# Planetary number
maxR = 0



data = pd.read_table('sullivan_table.txt', delim_whitespace=True, names=('RA', 'Dec', 'Rp', 'P', 'P insu', 'Rad v', 'Rs', 'Teff', 'Vmag', 'ICmag', 'Jmag', 'KSmag', 'Dmod', 'Dil p', 'devi flux', 'Sig-noi', 'NumPl'))

effTempSul = data['Teff']
sRad = data['Rs']
vMag = data['Vmag']
ICMag = data['ICmag']
jMag = data['Jmag']
KSMag = data['KSmag']
periodSul = data['P']
rPlanetSul = data['Rp']
RA = data['RA']
dec = data['Dec']

periodRange = 0.05		# ratio of period range
radiusRange = 0.05		# ratio of radius range


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
	effTempKep = [k.split(',')[26] for k in data]
	
	
pTESS = 0
rTESS = 0
count1 = 0
count2 = 0

outputFile = open('TESSData.csv', 'w')
for i in range(len(periodSul)): 		#
	for k in range(len(periodKep)):
		if rPlanetKep[k]:
			if effTempKep[k] > 4000:
				if periodSul[i] - periodRange*periodSul[i] <= float(periodKep[k]) <= periodSul[i] + periodRange*periodSul[i] and rPlanetSul[i] - radiusRange*rPlanetSul[i] <= float(rPlanetKep[k]) <= rPlanetSul[i] + radiusRange*rPlanetSul[i] and effTempSul[i] > 4000:
					ratioP = float(periodKep[k]) / periodSul[i]	
					ratioR = float(rPlanetKep[k]) / rPlanetSul[i]
					j = k
					count1 = 0
					count2 = 0
					countP = 0

						
					
					while kepID[j] == kepID[j+1]:
						if name[k][-2:] != '01' and count1 == 0:
							pNum = int(name[k][-2:])
							for l in range(1,pNum):
								pTESS = float(periodKep[k-l])*ratioP
								rTESS = float(rPlanetKep[k-l])*ratioR
								periodTESS.append(float(periodKep[k-l])/ratioP)
								rPlanetTESS.append(float(rPlanetKep[k-l])/ratioR)
								effTemp.append(float(effTempKep[k-l]))
								outputFile.write(name[k-l] + ',' +  kepID[k-l] + ',' + str(pTESS) + ',' + str(rTESS) + ',' + mStar[k-l] + ',' + numEpoch[k-l] + ','  + transitDur[k-l] + ',' + rStar[k-l] + ',' + str(RA[i]) + ',' + str(dec[i]) + ',' + effTempKep[k] + ',' + str(effTempSul[i]) + ',' + str(ICMag[i]) + '\n')
								countP = l
								

						pTESS = float(periodKep[k+count1])*ratioP
						rTESS = float(rPlanetKep[k+count1])*ratioR
						periodTESS.append(pTESS)
						rPlanetTESS.append(rTESS)
						effTemp.append(float(effTempKep[k]))
							
						if i + count1  >= len(RA):
							break
							
						j += 1		
						
						outputFile.write(name[k+count1] + ',' +  kepID[k+count1] + ',' + str(pTESS) + ',' + str(rTESS) + ',' + mStar[k+count1] + ',' + numEpoch[k+count1] + ','  + transitDur[k+count1] + ',' + rStar[k+count1] + ',' + str(RA[i]) + ',' + str(dec[i]) + ',' + effTempKep[k] + ',' + str(effTempSul[i]) + ',' + str(ICMag[i]) + '\n')
						if kepID[j] != kepID[j+1]:
							for n in range(0,int(count1) + 1 + int(countP)):
								numP.append(int(count1) + 1 + int(countP))

						if rTESS > maxR:
							#~ print rTESS, rPlanetKep[k+count1], pTESS, periodKep[k+count1], kepID[j]
							maxR = rTESS
							
							
						count1 += 1
						
						
							
		
					break
			elif effTempKep[k] <= 4000:
				if periodSul[i] - periodRange*periodSul[i] <= float(periodKep[k]) <= periodSul[i] + periodRange*periodSul[i] and rPlanetSul[i] - radiusRange*rPlanetSul[i] <= float(rPlanetKep[k]) <= rPlanetSul[i] + radiusRange*rPlanetSul[i] and effTempSul[i] < 4000:
					ratioP = float(periodKep[k]) / periodSul[i]	
					ratioR = float(rPlanetKep[k]) / rPlanetSul[i]
					j = k
					count1 = 0
					count2 = 0
					countP = 0
					
					while kepID[j] == kepID[j+1]:
						if name[k][-2:] != '01' and count1 == 0:
							pNum = int(name[k][-2:])
							for l in range(1,pNum):
								pTESS = float(periodKep[k-l])*ratioP
								rTESS = float(rPlanetKep[k-l])*ratioR	
								periodTESS.append(float(periodKep[k-l])/ratioP)
								rPlanetTESS.append(float(rPlanetKep[k-l])/ratioR)
								effTemp.append(float(effTempKep[k]))
								outputFile.write(name[k-l] + ' ' +  kepID[k-l] + ' ' + str(pTESS) + ' ' + str(rTESS) + ' ' + mStar[k-l] + ' ' + numEpoch[k-l] + ' '  + transitDur[k-l] + ' ' + rStar[k-l] + ' ' + str(RA[i+count1]) + ' ' + str(dec[i+count1]) + ' ' + effTempKep[k+count1] + ' ' + str(effTempSul[i+count1]) + ' ' + str(ICMag[i+count1]) + '\n')
								countP = l
								
						pTESS = float(periodKep[k+count1])*ratioP
						rTESS = float(rPlanetKep[k+count1])*ratioR
						periodTESS.append(pTESS)
						rPlanetTESS.append(rTESS)
						effTemp.append(float(effTempKep[k]))
							
						if i + count1  >= len(RA):
							break
							
						j += 1		
						
						outputFile.write(name[k+count1] + ' ' +  kepID[k+count1] + ' ' + str(pTESS) + ' ' + str(rTESS) + ' ' + mStar[k+count1] + ' ' + numEpoch[k+count1] + ' '  + transitDur[k+count1] + ' ' + rStar[k+count1] + ' ' + str(RA[i+count1]) + ' ' + str(dec[i+count1]) + ' ' + effTempKep[k+count1] + ' ' + str(effTempSul[i+count1]) + ' ' + str(ICMag[i+count1]) + '\n')
						if kepID[j] != kepID[j+1]:
							for n in range(0,int(count1) + 1 + int(countP)):
								numP.append(int(count1) + 1 + int(countP))

								
								
						if rTESS > maxR:
							#~ print rTESS, rPlanetKep[k+count1], pTESS, periodKep[k+count1]
							maxR = rTESS
						count1 += 1

						
						
		
					break
outputFile.close()
print len(rPlanetTESS), len(periodSul)					
periodTESS = np.log10(periodTESS)
rPlanetTESS = np.log10(rPlanetTESS)	
markers = ['.', 'o', '^', 's', 'p', 'h']
colors = ['black', 'blue', 'red', 'green', 'm', 'brown']



colorRange = np.linspace(0, 13, 13)


norm = mp.colors.Normalize(
    vmin=np.min(effTemp),
    vmax=np.max(effTemp))
    
c_m = mp.cm.cool
s_m = mp.cm.ScalarMappable(cmap=cm.plasma, norm=norm)
s_m.set_array([])

			
plt.scatter(periodTESS,rPlanetTESS,s=2, c = effTemp, cmap = cm.plasma )
plt.colorbar(s_m, label='Effective temperature of host star')
plt.ylim(-0.3,1.25)
plt.ylabel('log$_{10}$[Planet radius (R$_{\oplus}$)]', fontsize=12)
plt.xlabel('log$_{10}$[Period (days)', fontsize=12)
plt.savefig('plots/R_P-plot_effTemp1.pdf')
plt.clf()


for i in range(len(periodTESS)):
	if numP[i] != 1:
		markNum = (int(numP[i]) - 1)
		plt.scatter(periodTESS[i],rPlanetTESS[i], c = colors[markNum], marker=markers[markNum], edgecolor='black', label=numP, alpha=0.6)
plt.ylim(-0.3,1.25)


blue_circle = mlines.Line2D([], [], color='blue', marker='o', linestyle='None',
                          markersize=10, label='2')
red_triangle = mlines.Line2D([], [], color='red', marker='^', linestyle='None',
                          markersize=10, label='3')
green_square = mlines.Line2D([], [], color='green', marker='s', linestyle='None',
                          markersize=10, label='4')
m_pentagon = mlines.Line2D([], [], color='m', marker='p', linestyle='None',
                          markersize=10, label='5')
purple_hexagon = mlines.Line2D([], [], color='brown', marker='h', linestyle='None',
                          markersize=10, label='6')
                          
plt.legend(handles=[blue_circle, red_triangle, green_square, m_pentagon, purple_hexagon] ,loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=3, fancybox=True, shadow=True, title="Number of Planets")
plt.ylabel('log$_{10}$[Planet radius (R$_{\oplus}$)]', fontsize=12)
plt.xlabel('log$_{10}$[Period (days)]', fontsize=12)
plt.savefig('plots/R_P-plot_numP1.pdf')
plt.clf()



	
