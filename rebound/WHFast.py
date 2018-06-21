import rebound
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import math

start_time = time.time()

sim = rebound.Simulation()

sim.add(m=0.7740)

sim.add(m=2.20179294439e-05, a=0.12419, e=0.05157689260630252)

sim.add(m=2.14838440723e-05, a=0.125953, e=0.024280812017953325)

#~ sim.add(m=1.62353704347e-05, a=0.057657, e=0.03430694486301676)





sim.integrator = "whfast"
sim.dt = 1e-3
particles = sim.particles

torb = 2.*np.pi
Noutputs = 1000000
x1 = np.zeros(Noutputs)
y1 = np.zeros(Noutputs)
x2 = np.zeros(Noutputs)
y2 = np.zeros(Noutputs)
#~ x3 = np.zeros(Noutputs)
#~ y3 = np.zeros(Noutputs)
#~ x4 = np.zeros(Noutputs)
#~ y4 = np.zeros(Noutputs)
#~ x5 = np.zeros(Noutputs)
#~ y5 = np.zeros(Noutputs)

sim.move_to_com()

times = np.linspace(0, 900.*torb, Noutputs)
for i,t in enumerate(times):
    sim.integrate(t, exact_finish_time=0)
    x1[i] = particles[1].x
    y1[i] = particles[1].y
    x2[i] = particles[2].x
    y2[i] = particles[2].y
    #~ x3[i] = particles[3].x
    #~ y3[i] = particles[3].y
    #~ x4[i] = particles[4].x
    #~ y4[i] = particles[4].y
    #~ x5[i] = particles[5].x
    #~ y5[i] = particles[5].y


fig = plt.figure(figsize=(6,6))
ax = plt.subplot(111)
plt.scatter(x1, y1, marker='.', color='blue', s=1.2);
plt.scatter(x2, y2, marker='+', color='red', s=1.2);
#~ plt.scatter(x3, y3, marker='.', color='black', s=1.2);
#~ plt.scatter(x4, y4, marker='+', color='green', s=1.2);
#~ plt.scatter(x5, y5, marker='+', color='purple', s=1.2);
plt.scatter(0,0, marker='o', color='black', s=10);
plt.xlabel('Distance [AU]')
plt.ylabel('Distance [AU]')
plt.savefig('simResult.png')

execTime = open('execTime.txt', 'w')
execTime.write("--- %s seconds ---" % (time.time() - start_time) + '\n')
execTime.close()
