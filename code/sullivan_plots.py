import numpy as np
import matplotlib.pyplot as plt

effTemp = []		# effective temperature
sRad = [] 			# stellar radius
vMag = []			# V-band magnitude (all magnitudes are apparent magnitude)
ICMag = []			# I_C magnitude
jMag = []			# J-band magnitude
KSMag = []			# K_S magnitude


with open('sullivan_table.txt','r') as inputFile: 
	data = inputFile.readlines()[32:]
	sRad = [k.split(' ')[20] for k in data]
	effTemp = [k.split(' ')[22] for k in data]
	vMag = [k.split(' ')[24] for k in data]
	ICMag = [k.split(' ')[26] for k in data]
	jMag = [k.split(' ')[28] for k in data]
	KSMag = [k.split(' ')[31] for k in data]
	
print  sRad[0], effTemp[0],  vMag[0],  ICMag[0], jMag[0], KSMag[0]
print sRad[1], effTemp[1], vMag[1], ICMag[1], jMag[1], KSMag[1]

	
	
