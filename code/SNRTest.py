import math
import numpy as np


phCount = 3.958 * 10**11 * 10**(-0.4*float(12))
signalNoise = ((4*0.009158/1)**2)/(1/(math.sqrt(phCount) * 1/math.sqrt(6)))

print signalNoise

s = ((7.3*(1/math.sqrt(6)))/0.0013)**2
print s
mag = np.log10(s/(3.958*10**(11)))/-0.4
print mag
