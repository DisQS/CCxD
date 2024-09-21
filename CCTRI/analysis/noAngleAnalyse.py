import matplotlib.pyplot as plt
import numpy as np
import math
from pathlib import Path
from scipy.optimize import curve_fit
import os

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

steps = 36
size = 7
parallelamt = 0
analysesteps = 20
sym = 1
zrange = 25

thval = 100
phival = 0


cur_path = Path.cwd()
new_path = cur_path.parent /"Data" 
foldername = 'CCTRI-st' + str(steps)+'-si'+str(size + parallelamt)+ ('-sym' if sym==1 else '') + '-th' + str(thval) + '-ph' + str(phival)
if not (Path(cur_path/foldername).is_dir()):
    Path.mkdir(cur_path / foldername)
figz, axz = plt.subplots()
figt, axt = plt.subplots()
figth, axth = plt.subplots()
figg, axg = plt.subplots()

vars = ['$z$', '$t$', '$θ$', '$g$']
i = 0
for ax in axz, axt, axth, axg:
    ax.set_xlabel(vars[i])
    ax.set_ylabel('$P($' + vars[i] + '$)$')
    ax.set_title('$P($' + vars[i] + '$)$' + ' distribution for 10^' + str(size) + ' samples, ' + str(steps) + ' steps and ' + ('symmetrised' if sym==1 else 'not symmetrised'))
    ax.set_axisbelow(True)
    ax.grid(True)
    i+=1

axz.set_xlim(-10,10)
axz.set_ylim(0,0.275)

axg.set_xlim(0,1)
axg.set_ylim(0.5,2)

axt.set_xlim(0,1)
axt.set_ylim(0,3)

axth.set_xlim(0,1.571)
axth.set_ylim(0,1.2)






for i in range(0,analysesteps):
    z = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(thval) + "-" + str(phival)) / ("000") / ("dists" + str(i) + ".txt" )).readlines()


    thlength = int(z[0])
    tlength = int(z[1])
    glength = int(z[2])
    zlength = int(z[3])
    zstart = 4+int(thlength)+int(tlength)+int(glength)
    tstart = 4 + int(thlength)
    gstart = 4 + int(thlength) + int(tlength)

    zxdivide = zlength/(2*zrange)
    tydivide = (10**(size))/tlength
    zydivide = (10**size)/zxdivide

    thxvals = [th/tlength for th in range(thlength)]
    thyvals = [float(z[th].strip())/tydivide for th in range(4,4+thlength)]

    gyvals = [float(z[g].strip())/tydivide for g in range(gstart,gstart+glength)]

    txvals = [x/tlength for x in range(tlength)]
    tyvals = [float(z[t].strip())/tydivide for t in range(tstart,tstart+tlength)]
    zxvals = [(x/zxdivide)-zrange for x in range(zlength)]
    zyvals = [float(z[x].strip())/zydivide for x in range(zstart,zstart+zlength)]


    
    axz.scatter(zxvals,zyvals,s=5)
    axt.scatter(txvals,tyvals,s=5)
    axg.scatter(txvals,gyvals,s=5)
    axth.scatter(thxvals,thyvals,s=5)




figz.savefig(foldername + "/CCTRI-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-z.png")
figt.savefig(foldername + "/CCTRI-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-t.png")
figth.savefig(foldername + "/CCTRI-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-th.png")
figg.savefig(foldername + "/CCTRI-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-tg.png")

