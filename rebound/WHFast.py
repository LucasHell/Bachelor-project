import rebound
import numpy as np
import matplotlib.pyplot as plt
import time
import math

start_time = time.time()

sim = rebound.Simulation()

sim.add(m=1)

sim.add(m=1e-3, x=1., vy=1.)

sim.add(m=1e-2, a=4., e=0.1)


sim.integrator = "whfast"
sim.dt = 1e-3
particles = sim.particles

torb = 2.*np.pi
Noutputs = 5000
x1 = np.zeros(Noutputs)
y1 = np.zeros(Noutputs)
x2 = np.zeros(Noutputs)
y2 = np.zeros(Noutputs)

sim.move_to_com()

times = np.linspace(0, 100000, Noutputs)
for i,t in enumerate(times):
    sim.integrate(t, exact_finish_time=0)
    x1[i] = particles[1].x
    y1[i] = particles[1].y
    x2[i] = particles[2].x
    y2[i] = particles[2].y


fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
#~ ax.set_xlim([-1.5,1.5])
#~ ax.set_ylim([-1.5,1.5])
plt.scatter(x1, y1, marker='.', color='b', s=1.2);
plt.scatter(x2, y2, marker='+', color='r', s=1.2);
sim.status()
plt.savefig('simResult.png')

execTime = open('execTime.txt', 'w')
execTime.write("--- %s seconds ---" % (time.time() - start_time) + '\n')
execTime.close()
