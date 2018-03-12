import rebound
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()

sim.add(m=1)

sim.add(m=1e-3, x=1., vy=1.)

sim.add(m=1e-3, a=2., e=0.1)

sim.add(m=1e-3, a=2., e=0)

sim.integrator = "whfast"
sim.dt = 1e-3

sim.integrate(6.28318530717959, exact_finish_time=0)   # 6.28318530717959 is 2*pi

particles = sim.particles

torb = 2.*np.pi
Noutputs = 100
times = np.linspace(torb, 2.*torb, Noutputs)
x = np.zeros(Noutputs)
y = np.zeros(Noutputs)

for i,time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y
    
fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
plt.plot(x, y);
plt.show()


Noutputs = 1000
times = np.linspace(2.*torb, 20.*torb, Noutputs)
x = np.zeros(Noutputs)
y = np.zeros(Noutputs)
for i,time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
ax.set_xlim([-2,2])
ax.set_ylim([-2,2])
plt.plot(x, y);
plt.show()


sim.move_to_com()

times = np.linspace(20.*torb, 1000.*torb, Noutputs)
for i,time in enumerate(times):
    sim.integrate(time, exact_finish_time=0)
    x[i] = particles[1].x
    y[i] = particles[1].y

fig = plt.figure(figsize=(5,5))
ax = plt.subplot(111)
#~ ax.set_xlim([-1.5,1.5])
#~ ax.set_ylim([-1.5,1.5])
plt.scatter(x, y, marker='.', color='k', s=1.2);
plt.show()
