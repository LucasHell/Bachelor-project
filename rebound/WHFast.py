import rebound
import numpy as np
import matplotlib.pyplot as plt
import time
import math

start_time = time.time()

sim = rebound.Simulation()

sim.add(m=0.8560)

sim.add(m=1.378e-05, a=0.0542, e=0.0737)

sim.add(m=1.911e-05, a=0.087, e=0.0614)





sim.integrator = "whfast"
sim.dt = 1e-3
particles = sim.particles

torb = 2.*np.pi
Noutputs = 5000
x1 = np.zeros(Noutputs)
y1 = np.zeros(Noutputs)
x2 = np.zeros(Noutputs)
y2 = np.zeros(Noutputs)
x3 = np.zeros(Noutputs)
y3 = np.zeros(Noutputs)
#~ x4 = np.zeros(Noutputs)
#~ y4 = np.zeros(Noutputs)

sim.move_to_com()

times = np.linspace(0, 73446555.*torb, Noutputs)
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


fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
#~ ax.set_xlim([-1.5,1.5])
#~ ax.set_ylim([-1.5,1.5])
plt.scatter(x1, y1, marker='.', color='blue', s=1.2);
plt.scatter(x2, y2, marker='+', color='red', s=1.2);
#~ plt.scatter(x3, y3, marker='.', color='black', s=1.2);
#~ plt.scatter(x4, y4, marker='+', color='green', s=1.2);
#~ sim.status()
plt.savefig('simResult.png')

execTime = open('execTime.txt', 'w')
execTime.write("--- %s seconds ---" % (time.time() - start_time) + '\n')
execTime.close()
