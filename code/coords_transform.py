import numpy as np
import csv
import matplotlib.pyplot as plt
import math
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5


#~ RA = []
#~ dec = []
tilt = math.radians(23.439281)
#~ RA = np.linspace(0, 360, 360)
#~ dec = np.linspace(-90 , 90, 360)
beta = []
lamb = []
lamb2 = []

with open('TESSTime.csv', 'r') as inputFile:
	data = inputFile.readlines()
	RA = [k.split(',')[0] for k in data]
	dec = [k.split(',')[1] for k in data]

#~ RA = [50, 360, 280, 175]
#~ dec = [5, 0, 30, -90]

RA = map(float, RA)
dec = map(float, dec)

for n in range(len(RA)):
	RA[n] = math.radians(RA[n])
	dec[n] = math.radians(dec[n])
	
#~ for i in range(len(RA)):
	#~ for l in range(len(dec)):
		#~ print i, l
		#~ beta.append(math.asin(math.cos(tilt)*math.sin(dec[l]) - math.sin(RA[i])*math.cos(dec[l])*math.sin(tilt)))
		#~ print beta[i+l], RA[i], dec[l]
		#~ lamb.append(math.degrees((math.sin(tilt)*math.sin(dec[l]) + math.sin(RA[i])*math.cos(dec[l])*math.cos(tilt))/math.cos(beta[i+l])))
		#~ lamb2.append(math.degrees(math.acos((math.cos(RA[i])*math.cos(dec[l]))/math.cos(beta[i+l]))))
		#~ beta[i+l] = math.degrees(beta[i+l])

	
for i in range(len(RA)):
	beta.append(math.asin(math.cos(tilt)*math.sin(dec[i]) - math.sin(RA[i])*math.cos(dec[i])*math.sin(tilt)))
	lamb.append(math.degrees(math.asin((math.sin(tilt)*math.sin(dec[i]) + math.sin(RA[i])*math.cos(dec[i])*math.cos(tilt))/math.cos(beta[i]))))
	lamb2.append(math.degrees(math.acos((math.cos(RA[i])*math.cos(dec[i]))/math.cos(beta[i]))))
	beta[i] = math.degrees(beta[i])
	
for m in range(len(RA)):
	RA[m] = math.degrees(RA[m])
	dec[m] = math.degrees(dec[m])
	
#~ print lamb
#~ print lamb2
#~ print beta
print np.amax(dec), np.amin(dec)
print np.amax(RA), np.amin(RA)
print np.amax(lamb), np.amin(lamb)
print np.amax(lamb2), np.amin(lamb2)
c = SkyCoord(lon=lamb2*u.deg, lat=beta*u.deg, frame='heliocentrictrueecliptic')
ra_rad = c.lon.wrap_at(180 * u.deg).rad			
dec_rad = c.lat.rad

#~ print np.amax(lamb), np.amax(beta)

	
plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.scatter(ra_rad, dec_rad, s=7)
plt.savefig('plots/skymap_TESS_wrap_cutoff')
plt.clf()

plt.figure(figsize=(8,4.2))
plt.subplot(111, projection="aitoff")
plt.grid(True)
plt.title("Position of observed TESS objects", y=1.08)
plt.scatter(RA, dec, s=7)
plt.savefig('plots/skymap_RA_dec')
plt.clf()


print len(RA), len(dec), len(beta)
#~ plt.scatter(lamb, beta)
#~ plt.xlabel('Dec')
#~ plt.ylabel('Latitude')
#~ plt.savefig('plots/beta_RA.png')
#~ plt.clf()

#~ plt.scatter(RA, lamb)
#~ plt.xlabel('RA')
#~ plt.ylabel('Longitude')
#~ plt.savefig('plots/lamb_RA.png')
#~ plt.clf()

#~ plt.scatter(RA, lamb2)
#~ plt.xlabel('RA')
#~ plt.ylabel('Longitude')
#~ plt.savefig('plots/lamb2_RA.png')



