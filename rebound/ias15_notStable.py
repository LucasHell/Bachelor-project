import rebound
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations
import time as time2
import math

start_time = time2.time()

def setupSimulation():
	sim = rebound.Simulation()
	sim.integrator = "ias15" # IAS15 is the default integrator, so we don't need this line
	sim.add(m=0.7740)
	sim.add(m=2.20179294439e-05, a=0.12419, e=0.05157689260630252)
	sim.add(m=2.14838440723e-05, a=0.125953, e=0.024280812017953325)
	sim.move_to_com()
	return sim
    
def mergeParticles(sim):
    # Find two closest particles
    min_d2 = 1e9 # large number
    ps = sim.particles
    for i1, i2 in combinations(range(sim.N),2): # get all pairs of indices
        dp = ps[i1] - ps[i2]   # Calculates the coponentwise difference between particles
        d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
        if d2<min_d2:
            min_d2 = d2
            col_i1 = i1
            col_i2 = i2

    cp1 = ps[col_i1]
    cp2 = ps[col_i2]
    # Merge two closest particles

    sum_mass = cp1.m + cp2.m
    mergedPlanet = (cp1*cp1.m + cp2*cp2.m)/sum_mass
    mergedPlanet.m  = sum_mass
    sim.remove(index=col_i2)
    sim.remove(index=col_i1)
    sim.add(mergedPlanet, assignHash=True)


sim = setupSimulation() # Resets everything
sim.exit_min_distance = 0.005
Noutputs = 100000
times = np.linspace(0,900.*2.*np.pi,Noutputs)
distances = np.zeros(Noutputs)
x1 = np.zeros(Noutputs)
y1 = np.zeros(Noutputs)
x2 = np.zeros(Noutputs)
y2 = np.zeros(Noutputs)
#~ x3 = np.zeros(Noutputs)
#~ y3 = np.zeros(Noutputs)
ps = sim.particles # ps is now an array of pointers. It will update as the simulation runs.
try:
    for i,time in enumerate(times):
		sim.integrate(time)
		dp = ps[1] - ps[2]   # Calculates the coponentwise difference between particles
		distances[i] = np.sqrt(dp.x*dp.x+dp.y*dp.y+dp.z*dp.z)
		x1[i] = ps[1].x
		y1[i] = ps[1].y
		x2[i] = ps[2].x
		y2[i] = ps[2].y
		#~ x3[i] = ps[3].x
		#~ y3[i] = ps[3].y
except rebound.Encounter as error:
    print(error)
    
    
print("Number of particles at the beginning of the simulation: %d."%sim.N)
for i,time in enumerate(times):
    try:
        sim.integrate(time)
    except rebound.Encounter as error:
        print(error)
        mergeParticles(sim)
print("Number of particles at the end of the simulation: %d."%sim.N)

#~ fig = plt.figure(figsize=(10,5))
#~ ax = plt.subplot(111)
#~ ax.set_xlabel("time [orbits]")
#~ ax.set_xlim([0,sim.t/(2.*np.pi)])
#~ ax.set_ylabel("distance")
#~ plt.plot(times/(2.*np.pi), distances);
#~ plt.plot([0.0,12],[0.2,0.2]); # Plot our close encounter criteria;
#~ plt.clf()

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
#~ ax.set_xlim([-1.5,1.5])
#~ ax.set_ylim([-1.5,1.5])
plt.scatter(x1, y1, marker='.', color='blue', s=1.2);
plt.scatter(x2, y2, marker='+', color='red', s=1.2);
#~ plt.scatter(x3, y3, marker='.', color='black', s=1.2);
#~ plt.scatter(x4, y4, marker='+', color='green', s=1.2);
#~ plt.scatter(x5, y5, marker='+', color='purple', s=1.2);
plt.scatter(0,0, marker='o', color='black', s=10);
plt.savefig('simResult.png')

execTime = open('execTime.txt', 'w')
execTime.write("--- %s seconds ---" % (time2.time() - start_time) + '\n')
execTime.close()


