import matplotlib.pyplot as plt
import csv

epoch1 = []
transitTime1 = []
rSky1 = []
vSky1 = []
epoch2 = []
transitTime2 = []
rSky2 = []
vSky2 = []

with open('Times','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
		if row[0] == 0:
			print(test)
			epoch1.append(long(row[1]))
			transitTime1.append(long(row[2]))
			rSky1.append(long(row[3]))
			vSky1.append(long(row[4]))
	#~ else:
		#~ epoch2.append(float(row[1]))
		#~ transitTime2.append(float(row[2]))
		#~ rSky2.append(float(row[3]))
		#~ vSky2.append(float(row[4]))
		
f = open('Times')
yourList = f.readlines()

#print(yourList)
print(plots)
		
			
#plt.plot(epoch1,transitTime1, label='Mass of exoplanet')
#plt.xlabel('Planet Radius $R_E$')
#plt.ylabel('Planet Mass $M_E$')
#plt.title('Mass-radius relation')
#plt.legend()
#plt.show()

