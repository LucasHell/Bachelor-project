import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib  import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import matplotlib.lines as mlines
import math

periodFrac = []
ampl = []


with open('AmplPeriod.csv','r') as inputFile: 
	data = inputFile.readlines()
	periodFrac = [k.split(',')[0] for k in data]
	ampl = [k.split(',')[1] for k in data]
	
	
plt.plot(periodFrac,ampl)
plt.show()

