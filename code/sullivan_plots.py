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

data = pd.read_table('sullivan_table.txt', delim_whitespace=True, names=('RA', 'Dec', 'Rp', 'P', 'P insu', 'Rad v', 'Rs', 'Teff', 'Vmag', 'ICmag', 'Jmag', 'KSmag', 'Dmod', 'Dil p', 'devi flux', 'Sig-noi', 'NumPl'))

effTemp = data['RA']
sRad = data['Rs']
vMag = data['Vmag']
ICMag = data['ICmag']
jMag = data['Jmag']
KSMag = data['KSmag']

plt.hist(effTemp,bins=50)
plt.title("Histogram of the effective temperatures of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/effTemp.png')
plt.clf()

plt.hist(sRad,bins=50)
plt.title("Histogram of the star radii of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/sRad.png')
plt.clf()

plt.hist(vMag,bins=50)
plt.title("Histogram of the V-band magnitude of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/vMag.png')
plt.clf()

plt.hist(ICMag,bins=50)
plt.title("Histogram of the I_C_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/ICMag.png')
plt.clf()

plt.hist(jMag,bins=50)
plt.title("Histogram of the J-band magnitude of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/jMag.png')
plt.clf()

plt.hist(KSMag,bins=50)
plt.title("Histogram of the K_s_ band magnitude of objects\nfrom the Sullivan catalogue")
plt.savefig('./plots/histo/KSMag.png')



	
	
