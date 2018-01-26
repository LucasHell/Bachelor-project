import matplotlib.pyplot as plt
import numpy as np



# Data for plotting
R = np.linspace(0.0, 4.0, 40)
mass = []

def frange(start,stop, step=0.1):
    while start < stop:
        yield start
        start +=step

for i in frange(0,4,0.1):
	if i < 1.5:
		M = 0.440*(i**3) + 0.614*(i**4);
		print(M, i)
		mass.append(M)
	else: 
		M = 2.69*(i**0.93)
		print(M, i)
		mass.append(M)


#Note that using plt.subplots below is equivalent to using
#fig = plt.figure and then ax = fig.add_subplot(111)
fig, ax = plt.subplots()
ax.plot(R, mass)

ax.set(xlabel='Planet Radius ($R_E$)', ylabel='Planet Mass ($M_E$)')

fig.savefig("massRadius.png")
plt.show()
