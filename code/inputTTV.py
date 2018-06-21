import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import pathlib2
import pandas as pd
from fractions import Fraction  
import math

G = 0.000295994511
inputFile = []
data = []
name = [] 		# KOI name of planet
kepID = [] 		# kepler ID of planet
mStar = [] 		# mass of host star
mPlanet = [] 	# mass of planet
period = []		# period of planet
rPlanet = []	# radius of planet
numEpoch = []	# number of epochs
eccentricity = 0	# eccentricity of planet
inclination = 90	# icnlination of planet
lNode = 0			# long node of planet
argument = 0		# argument of planet
meanAnom = []		# mean anomaly of planet
count = 0			# counter for number of planets
kepMag = []			# Kepler magnitude of star
transitDur = []		# Duration of transit
rStar = []			# Radius of star
RA = []				# Right ascension of system
dec = []			# Declination of system
ofirName = []
ofirID = []
ofirPeriod = []
ofirRPlanet = []
ofirMStar = []
ofirEpoch = []
ofirDur = []
ofirRStar = []
ofirRA = []
ofirDec = []
ofirEffTemp = []
ofirKepMag = []
sumMass = []
sysmPlanet = []
sysPeriod = []
sysEccentricty = []
sysInclination = []
syslNode = []
sysArg = []
sysMeanAnom = []
sysRPlanet = []
sysCountList = []




# read in data from csv file to respective arrays
#~ with open('nasaDATA.csv','r') as inputFile:
	#~ data = inputFile.readlines()[46:]
	#~ name = [k.split(',')[2] for k in data]
	#~ kepID = [k.split(',')[1] for k in data]
	#~ period = [k.split(',')[5] for k in data]
	#~ rPlanet = [k.split(',')[20] for k in data]
	#~ mStar = [k.split(',')[32] for k in data]
	#~ numEpoch = [k.split(',')[8] for k in data]
	#~ kepMag = [k.split(',')[37] for k in data]
	#~ transitDur = [k.split(',')[14] for k in data]
	#~ rStar = [k.split(',')[29] for k in data]
	
	
# read in data from csv file to respective arrays
with open('TESSData.csv','r') as inputFile: 
	data = inputFile.readlines()
	name = [k.split(',')[0] for k in data]
	kepID = [k.split(',')[1] for k in data]
	period = [k.split(',')[2] for k in data]
	rPlanet = [k.split(',')[3] for k in data]
	mStar = [k.split(',')[4] for k in data]
	numEpoch = [k.split(',')[5] for k in data]
	transitDur = [k.split(',')[6] for k in data]
	rStar = [k.split(',')[7] for k in data]
	RA = [k.split(',')[8] for k in data]
	dec = [k.split(',')[9] for k in data]
	effTemp = [k.split(',')[10] for k in data]
	kepMag = [k.split(',')[12] for k in data]
	
#~ with open('nasa_ofir.csv','r') as inputFile: 
	#~ data = inputFile.readlines()[46	:]
	#~ name = [k.split(',')[2] for k in data]
	#~ kepID = [k.split(',')[1] for k in data]
	#~ period = [k.split(',')[5] for k in data]
	#~ rPlanet = [k.split(',')[20] for k in data]
	#~ mStar = [k.split(',')[32] for k in data]
	#~ numEpoch = [k.split(',')[8] for k in data]
	#~ transitDur = [k.split(',')[15] for k in data]
	#~ rStar = [k.split(',')[29] for k in data]
	#~ RA = [k.split(',')[35] for k in data]
	#~ dec = [k.split(',')[36] for k in data]
	#~ effTemp = [k.split(',')[24] for k in data]
	#~ kepMag = [k.split(',')[37] for k in data]
	
	
#~ data = pd.read_table('ofir_table.txt', sep=';', skiprows=34, names=('KOI_num', 'newDetFlag', 'TTVfre', 'TTV+uncer', 'TTV-uncer', 'TTV_per', 'Delta_chi', 'chi_area', 'chi_single', 'chi_RMS', 'cho_correl', 'TTV_amp', 'TTV_amp+_uncer', 'TTV_amp-_uncer', 'TTV_ref', 'TTV_ref+_uncer', 'TTV_ref-_uncer', 'cofid', 'STD_error', '20', '21', '22', '23', '24'))
#~ name_ofir = data['KOI_num']
#~ print name_ofir[8], name_ofir[9]
#~ noCount = 0

#~ for n in range(len(name_ofir)):
	#~ for m in range(len(name)):
		#~ if str(name_ofir[n]) == str(name[m][-len(str(name_ofir[n])):]):
			#~ if str(name[m][-2:]) != '01' and kepID[m] != kepID[m-1]:
				#~ for l in range(0,int(name[m][-2:])):
					#~ ofirName.append(name[m-l])
					#~ ofirID.append(kepID[m-l])
					#~ ofirPeriod.append(period[m-l])
					#~ ofirRPlanet.append(rPlanet[m-l])
					#~ ofirMStar.append(mStar[m-l])
					#~ ofirEpoch.append(numEpoch[m-l])
					#~ ofirDur.append(transitDur[m-l])
					#~ ofirRStar.append(rStar[m-l])
					#~ ofirRA.append(RA[m-l])
					#~ ofirDec.append(dec[m-l])
					#~ ofirEffTemp.append(effTemp[m-l])
					#~ ofirKepMag.append(kepMag[m-l])
				#~ break
					
			#~ else:

				#~ ofirName.append(name[m])
				#~ ofirID.append(kepID[m])
				#~ ofirPeriod.append(period[m])
				#~ ofirRPlanet.append(rPlanet[m])
				#~ ofirMStar.append(mStar[m])
				#~ ofirEpoch.append(numEpoch[m])
				#~ ofirDur.append(transitDur[m])
				#~ ofirRStar.append(rStar[m])
				#~ ofirRA.append(RA[m])
				#~ ofirDec.append(dec[m])
				#~ ofirEffTemp.append(effTemp[m])
				#~ ofirKepMag.append(kepMag[m])
				#~ break
		#~ if m == len(name)-1:
			#~ print "no match: ", name_ofir[n]
			#~ noCount += 1


#~ print len(name), len(ofirName), len(name_ofir), noCount

#~ name = ofirName
#~ kepID = ofirID
#~ period = ofirPeriod
#~ rPlanet = ofirRPlanet
#~ mStar = ofirMStar
#~ numEpoch = ofirEpoch
#~ transitDur = ofirDur
#~ rStar = ofirRStar
#~ RA = ofirRA
#~ dec = ofirDec
#~ effTemp = ofirEffTemp
#~ kepMag = ofirKepMag





#Convert planet radius to mass


for i in range(len(rPlanet)):
	if not rPlanet[i]:
		rPlanet[i] = 0
	else:
		rPlanet[i] = float(rPlanet[i]) 										# convert strings in rPlanet to floats
	
	if rPlanet[i] < 1.5: 													# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(0.440*(rPlanet[i]**3) + 0.614*(rPlanet[i]**4))
	else:																	# calculate mass of planet for planets with mass below 1.5 earth masses
		mPlanet.append(2.69*(rPlanet[i]**0.93))

	mPlanet[i] = mPlanet[i]*0.000002988 		#convert planet mass from earth masses to solar masses
#~ print np.amax(rPlanet)	
pathlib2.Path('./input/').mkdir(parents=True, exist_ok=True)				# create directory input if it does not exist

for i in range(len(dec)):
	if float(dec[i]) > 0:
		dec[i] = float(dec[i]) * -1
		dec[i] = str(dec[i])
		RA[i] = float(RA[i]) * -1
		RA[i] = str(RA[i])

# open output files to write data for TTVFast
systemCount = 0	
outputFile = open('input/%s.in' % systemCount, 'w')
outNumP = open('numberPlanets.csv', 'w')
outErrorFile = open('timingErrors.csv', 'w')
outTESSTime = open('TESSTime.csv', 'w')
outStabSim = open('stabSim.csv', 'w')
outRAdec = open('RA_dec_sys.csv', 'w')
outRPlanet = open('rPlanet.csv', 'w')
outputFile.write(repr(G) + '\n' + mStar[0] + '\n')




totalPlanets = 0
semiMajor = 0
rHill = 0
sumMass = []
semiMajorList = []
semiMajorNe = 0
dif = []
posMin = 0

#~ # write to output files in the format required for TTVFast
for i in range(len(mPlanet)): 			
	argument = np.random.uniform(0,360)
	meanAnom.append(90 - 360 * (float(numEpoch[i]) / float(period[i])) - argument)		# calculate mean anomaly of planet with reference to transit of first planet
	while meanAnom[i] > 360:
		meanAnom[i] = meanAnom[i] - 360 									# if angle of mean anomaly is above 360 degrees, subtract 360 until it is in the range 0 to 360 degrees
		if meanAnom[i] == 360:
			meanAnom[i] = 0
				
	while meanAnom[i] < 0:													# if angle of mean anomaly is below 360 degrees, add 360 until it is in the range 0 to 360
		meanAnom[i] = meanAnom[i] + 360
		if meanAnom[i] == 360:
			meanAnom[i] = 0
			
	if not rPlanet[i]:
		continue

	# write data for TTVFast
	#~ outputFile.write(repr(mPlanet[i]) + '\n' + period[i] + ' ' + str(np.random.rayleigh(0.03)) + ' ' + repr(inclination) + ' ' + repr(lNode) + '  ' + str(np.random.uniform(0,360)) + ' ' + repr(meanAnom[i]) + '\n') 
	sysmPlanet.append(mPlanet[i])
	sysPeriod.append(float(period[i]))
	sysEccentricty.append(np.random.rayleigh(0.03))
	sysInclination.append(inclination)
	syslNode.append(lNode)
	sysArg.append(argument)
	sysMeanAnom.append(meanAnom[i])
	sysRPlanet.append(rPlanet[i])
	count += 1			# counter for number of planets for each system
	
	
	if i == len(mPlanet)-1: 			# for last element
		if count == 1:
			os.remove('input/%s.in' % systemCount)
			outputFile.close()	
	
			
			
		elif count != 1:
			outNumP.write(repr(count) + '\n')
			S = 3.958 * 10**11 * 10**(-0.4*float(kepMag[i]))			# Kep: 7.8 * 10**8 *10**(-0.4*float(kepMag[i]))		TESS: 3.96 * 10**13 * 10**(-0.4*float(kepMag[i]))
			
			errorTiming = ((S * float(transitDur[i]))**Fraction('-1/2') * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**Fraction('-3/2') * float(transitDur[i]))		# error in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			

			for n in range(0,count):
				outTESSTime.write(RA[i] + ',' + dec[i] + '\n')
				semiMajor = ((float(period[i-n])**2 * G * float(mStar[i]))/(4 * math.pi**2))**Fraction('1/3')
				semiMajorList.append(semiMajor)
				sumMass.append(mPlanet[i-n])
				posMin = sysPeriod.index(np.amin(sysPeriod))
				outputFile.write(str(sysmPlanet[posMin]) + '\n' + str(sysPeriod[posMin]) + ' ' + repr(sysEccentricty[posMin]) + ' ' + repr(sysInclination[posMin]) + ' ' + repr(syslNode[posMin]) + '  ' + repr(sysArg[posMin]) + ' ' + repr(sysMeanAnom[posMin]) + '\n')
				sysPeriod[posMin] = 1000000
				outRPlanet.write(str(sysRPlanet[posMin]) + '\n')
			semiMajorList = sorted(semiMajorList, key=float, reverse=False)
			
			rHill = (semiMajorList[0] + semiMajorList[1])/2 * ((sumMass[0] + sumMass[1])/(3*float(mStar[i]))**Fraction('1/3'))
			dif.append((semiMajorList[1] - semiMajorList[0])/rHill)
			sysCountList.append(systemCount)
			outStabSim.write(str(rHill) + ',' + str(systemCount) +  '\n')
			
			totalPlanets += count	
			outputFile.close()
			outRAdec.write(RA[i] + ',' + dec[i] + '\n')
			break

			

	elif kepID[i] != kepID[i+1]:		# if ID is not the same as ID of next planet a new file is created for new system
		if count == 1:					# if number of planets is one the file is removed
			os.remove('input/%s.in' % systemCount)
			outputFile.close()	


			# begin on new file
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			count = 0
			semiMajorList = []
			sysmPlanet = []
			sysPeriod = []
			sysEccentricty = []
			sysInclination = []
			syslNode = []
			sysArg = []
			sysMeanAnom = []
			sysRPlanet = []
			
			
		elif count != 1:				# if number of planets is not 1 the rest of the data is saved and a new file is created
			outNumP.write(repr(count) + '\n')
			S = 3.958 * 10**11 * 10**(-0.4*float(kepMag[i]))			# Kep: 7.8 * 10**8 *10**(-0.4*float(kepMag[i]))		TESS: 3.96 * 10**13 * 10**(-0.4*float(kepMag[i]))
			errorTiming = ((S * float(transitDur[i]))**Fraction('-1/2') * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**Fraction('-3/2') * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			for n in range(0,count):
				outTESSTime.write(RA[i] + ',' + dec[i] + '\n')
				semiMajor = ((float(period[i-n])**2 * G * float(mStar[i]))/(4 * math.pi**2))**Fraction('1/3')
				semiMajorList.append(semiMajor)
				sumMass.append(mPlanet[i-n])	
				posMin = sysPeriod.index(np.amin(sysPeriod))
				outputFile.write(str(sysmPlanet[posMin]) + '\n' + str(sysPeriod[posMin]) + ' ' + repr(sysEccentricty[posMin]) + ' ' + repr(sysInclination[posMin]) + ' ' + repr(syslNode[posMin]) + '  ' + repr(sysArg[posMin]) + ' ' + repr(sysMeanAnom[posMin]) + '\n')
				sysPeriod[posMin] = 1000000
				outRPlanet.write(str(sysRPlanet[posMin]) + '\n')
			semiMajorList = sorted(semiMajorList, key=float, reverse=False)
			

			rHill = (semiMajorList[0] + semiMajorList[1])/2 * ((sumMass[0] + sumMass[1])/(3*float(mStar[i]))**Fraction('1/3'))
			dif.append((semiMajorList[1] - semiMajorList[0])/rHill)
			sysCountList.append(systemCount)
			outStabSim.write(str(rHill) + ',' + str(systemCount) +  '\n')
			
			totalPlanets += count	
			count = 0
			sumMass = []
			semiMajor = 0
			rHill = 0
			systemCount += 1
			semiMajorList = []
			sysmPlanet = []
			sysPeriod = []
			sysEccentricty = []
			sysInclination = []
			syslNode = []
			sysArg = []
			sysMeanAnom = []
			sysRPlanet = []
			
			outputFile.close()
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')	
			outRAdec.write(RA[i] + ',' + dec[i] + '\n')
			
			
	elif kepID[i] == kepID[i+1] and int(str(name[i])[-2:]) > int(str(name[i+1])[-2:]):
		if count == 1:					# if number of planets is 1 the file is removed
			os.remove('input/%s.in' % systemCount)
			outputFile.close()	
			systemCount += 1
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')
			count = 0

			
			
		elif count != 1:				# if number of planets is not 1 the rest of the data is saved and a new file is created
			outNumP.write(repr(count) + '\n')
			S = 3.958 * 10**11 * 10**(-0.4*float(kepMag[i]))			# Kep: 7.8 * 10**8 *10**(-0.4*float(kepMag[i]))		TESS: 3.96 * 10**13 * 10**(-0.4*float(kepMag[i]))
			errorTiming = ((S * float(transitDur[i]))**Fraction('-1/2') * ((float(rPlanet[i])*0.009158)/float(rStar[i]))**Fraction('-3/2') * float(transitDur[i]))		# timing precision in hours
			outErrorFile.write(repr(errorTiming*60) + '\n') 		# write errors to file in minutes
			
			for n in range(0,count):
				outTESSTime.write(RA[i] + ',' + dec[i] + '\n')
				semiMajor = ((float(period[i-n])**2 * G * float(mStar[i]))/(4 * math.pi**2))**Fraction('1/3')
				semiMajorList.append(semiMajor)
				sumMass.append(mPlanet[i-n])
				posMin = sysPeriod.index(np.amin(sysPeriod))
				outputFile.write(str(sysmPlanet[posMin]) + '\n' + str(sysPeriod[posMin]) + ' ' + repr(sysEccentricty[posMin]) + ' ' + repr(sysInclination[posMin]) + ' ' + repr(syslNode[posMin]) + '  ' + repr(sysArg[posMin]) + ' ' + repr(sysMeanAnom[posMin]) + '\n')
				sysPeriod[posMin] = 1000000
				outRPlanet.write(str(sysRPlanet[posMin]) + '\n')
			semiMajorList = sorted(semiMajorList, key=float, reverse=False)

			rHill = (semiMajorList[0] + semiMajorList[1])/2 * ((sumMass[0] + sumMass[1])/(3*float(mStar[i]))**Fraction('1/3'))
			dif.append((semiMajorList[1] - semiMajorList[0])/rHill)
			sysCountList.append(systemCount)
			outStabSim.write(str(rHill) + ',' + str(systemCount) +  '\n')
			
			totalPlanets += count	
			count = 0
			sumMass = []
			semiMajor = 0
			rHill = 0
			systemCount += 1
			semiMajorList = []
			sysmPlanet = []
			sysPeriod = []
			sysEccentricty = []
			sysInclination = []
			syslNode = []
			sysArg = []
			sysMeanAnom = []
			sysRPlanet = []
			
			outputFile.close()
			outputFile = open('input/%s.in' % systemCount, 'w')
			outputFile.write(repr(G) + '\n' + mStar[i+1] + '\n')	
			outRAdec.write(RA[i] + ',' + dec[i] + '\n')
			
print totalPlanets			
outDif = open('dif_table.csv', 'w')
for l in range(len(dif)):
	outDif.write(str(dif[l]) + ',' + str(sysCountList[l]) + ',' + '\n')
	
outDif.close()	
		
dif = sorted(dif, key=float, reverse=True)

outDif = open('dif_table_sort.csv', 'w')
for l in range(len(dif)):
	outDif.write(str(dif[l]) + ',' + str(sysCountList[l]) + ',' + '\n')
	
outNumP.close()
outErrorFile.close()
outTESSTime.close()
outStabSim.close()



	




