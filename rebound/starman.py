import rebound
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()
#~ sim.add(["Sun","Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune"],date="2018-02-10 00:00")
#~ sim.save("ss.bin")

sim = rebound.Simulation.from_file("ss.bin")
sim.add(primary=sim.particles[0],
        M=3.68007 *np.pi/180.,
        a=1.34126487,
        omega = 177.28664 *np.pi/180.,
        Omega = 317.45885 *np.pi/180.,
        e = 0.2648281,
        inc = 1.09424 *np.pi/180.)
        
tesla = sim.particles[-1]
earth = sim.particles[3]
r=np.linalg.norm(np.array(tesla.xyz) - np.array(earth.xyz))
v=np.linalg.norm(np.array(tesla.vxyz) - np.array(earth.vxyz))
energy = 0.5*v*v-earth.m/r
c3 = 2.*energy*887.40652 # from units where G=1, length=1AU to km and s
print("c3 = %f (km^2/s^2)" % c3)

rebound.OrbitPlot(sim,lim=1.8,slices=True,color=True);

sim.dt = sim.particles[1].P/60. # small fraction of Mercury's period
sim.integrator="mercurius"
N = 1000
times = np.linspace(0.,2.*np.pi*1e5,N)
a = np.zeros(N)
e = np.zeros(N)
for i,t in enumerate(times):
    sim.integrate(t,exact_finish_time=0)
    orbit = sim.particles[-1].calculate_orbit(primary=sim.particles[0])
    a[i] = orbit.a
    e[i] = orbit.e
    


fig = plt.figure(figsize=(9,7))
ax = plt.subplot(211)
ax.set_xlim([0,np.max(times)/2./np.pi])
ax.set_xlabel("time [yrs]")
ax.set_ylabel("semi-major axis [AU]")
plt.plot(times/2./np.pi,a)
ax = plt.subplot(212)
ax.set_xlim([0,np.max(times)/2./np.pi])
ax.set_xlabel("time [yrs]")
ax.set_ylabel("eccentricity")
plt.plot(times/2./np.pi,e);


