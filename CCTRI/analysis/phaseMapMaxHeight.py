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

steps = 40
size = 6
sym = 0
zrange = 25


cur_path = Path.cwd()
new_path = cur_path.parent / ("TRIRG-"  + str(size) + "-" + str(steps) + "-" + str(sym)) / "Data"



thmin = 0
thmax = 50
thstep = 1

psimin = 0
psimax = 50
psistep = 1

zmaxheight = [[0 for x in range(math.ceil((thmax-thmin)/thstep))] for y in range(math.ceil((psimax-psimin)/psistep))]
stddevs = [[0 for x in range(math.ceil((thmax-thmin)/thstep))] for y in range(math.ceil((psimax-psimin)/psistep))]

for th in range(thmin,thmax,thstep):

    for psi in range(psimin,psimax,psistep):

        i = math.floor((th-thmin)/thstep)
        j = math.floor((psi-psimin)/psistep)
        try:
            z = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(th) + "-" + str(psi)) / "dists" + str(steps) + ".txt" ).readlines()
        except FileNotFoundError:
            z = [0]
        if(len(z)>1):
            thlength = z[0]
            tlength = z[1]
            glength = z[2]
            zlength = z[3]
            zstart = 4+thlength+tlength+glength
            xdivide = zlength/(2*zrange)
            xvals = [(x/xdivide)-zrange for x in range(len(zlength))]
            yvals = [int(z[x].strip()) for x in range(zstart,zstart+zlength)]
            mean = sum(xvals*yvals)/zlength
            sigma = sum(yvals*(xvals-mean)**2)/zlength

            popt,pcov = curve_fit(gaus,xvals,yvals,p0=[1,mean,sigma])
            zmaxheight[i][j] = popt[2]
            stddevs[i][i] = popt[3]
        else:
            xvals = [0]
            yvals = [0]
            mean=0
            sigma=0
            zmaxheight[i][j] = 0

plt.figure("phasemap")
plt.imshow(zmaxheight,cmap="hot",interpolation="nearest")
plt.savefig("phasemap-" + str(size) + "-" + str(steps) + "-" + str(sym) + ".png")

plt.figure("phasemapSTDDEVS")
plt.imshow(stddevs,cmap="hot",interpolation="nearest")
plt.savefig("phasemapSTDDEVS-" + str(size) + "-" + str(steps) + "-" + str(sym) + ".png")

save = open("phasemap-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-ZMAX.txt","w")
save.write(zmaxheight)
save.close()

save = open("phasemap-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-STDDEVS.txt","w")
save.write(stddevs)
save.close()
        
