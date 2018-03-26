import rebound
import numpy as np
import matplotlib.pyplot as plt
import time
import math

start_time = time.time()

sim = rebound.Simulation()

sim.add(m=1)

sim.add(m=1e-3, x=1., vy=1.)

sim.add(m=1e-3, a=2., e=0.1)


sim.integrator = "whfast"
sim.dt = 1e-3
particles = sim.particles

torb = 2.*np.pi
Noutputs = 1000
x = np.zeros(Noutputs)
y = np.zeros(Noutputs)

sim.move_to_com()

times = np.linspace(0, 10000, Noutputs)
for i,t in enumerate(times):
    sim.integrate(t, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
#~ ax.set_xlim([-1.5,1.5])
#~ ax.set_ylim([-1.5,1.5])
plt.scatter(x, y, marker='.', color='k', s=1.2);
sim.status()
plt.savefig('simResult.png')

execTime = open('execTime.txt', 'w')
execTime.write("--- %s seconds ---" % (time.time() - start_time) + '\n')
execTime.close()
