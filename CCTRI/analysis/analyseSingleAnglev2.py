import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
import math
from pathlib import Path
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

steps = 12
size = 6
sym = 0
zrange = 25

thval = 0
phival = 0


cur_path = Path.cwd()
new_path = cur_path.parent /"Data"


z = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(thval) + "-" + str(phival)) / ("dists" + str(steps) + ".txt" )).readlines()


thlength = int(z[0])
tlength = int(z[1])
glength = int(z[2])
zlength = int(z[3])
zstart = 4+int(thlength)+int(tlength)+int(glength)
xdivide = zlength/(2*zrange)
xvals = [(x/xdivide)-zrange for x in range(zlength)]
yvals = [float(z[x].strip()) for x in range(zstart,zstart+zlength)]
mean = sum(xvals*yvals)/zlength
sigma = sum(yvals*(xvals-mean)**2)/zlength

popt,pcov = curve_fit(gaus,xvals,yvals,p0=[1,mean,sigma])

plt.figure("z")
plt.scatter(xvals,yvals,s=5)
plt.savefig("testz.png")
