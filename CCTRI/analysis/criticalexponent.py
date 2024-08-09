
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
import math
from pathlib import Path
import pylab as plb
from scipy.optimize import curve_fit

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


steps = 33
size = 7
analysesteps = 8
sym = 0 
zrange = 25

cur_path = Path.cwd()
new_path = cur_path.parent / "Data"

# array of strings to iterate over
offsetvalues = []


rayplot = [[0 for i in range(len(offsetvalues))] for j in range(analysesteps)]

critexp = [0 for i in range(analysesteps)];

for i in range(len(offsetvalues)):
    z = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(thval) + "-" + str(phival)) / ("dists" + str(i) + ".txt" )).readlines()
    z
