import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
import math
from pathlib import Path
import pylab as plb
from scipy.optimize import curve_fit
from numpy import random

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x,m,c):
    return (m * x) + c

steps = 20
size = 6 
analysesteps = 12
sym = 1 
zrange = 25
topxpercent = 10
subSampleAmtForGaussApprox = 10



cur_path = Path.cwd()
new_path = cur_path.parent /"Data"

pushval=["010","030","040","050"]
z = [[None for y in range(analysesteps)] for x in range(len(pushval))]
avgmean = [[0.0 for y in range(len(pushval))] for x in range(analysesteps)]
minmean = [[50.0 for y in range(len(pushval))] for x in range(analysesteps)]
maxmean = [[-50.0 for y in range(len(pushval))] for x in range(analysesteps)]


for i in range(0,analysesteps):
    for j in range(0,len(pushval)):
        zopen = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(sym) ) / pushval[j] / ("dists" + str(i) + ".txt" )).readlines()

        #slog through the data and find out where the z values actually start
        thlength = int(zopen[0])
        tlength = int(zopen[1])
        glength = int(zopen[2])
        zlength = int(zopen[3])
        zstart = 4+int(thlength)+int(tlength)+int(glength)
        print(str(thlength)+ ' ' +  str(tlength) + ' ' +  str(glength)  + ' ' + str(zlength))
        print(zstart)
        print(zlength)
        #constant to rescale and normalise data
        xdivide = zlength/(2*zrange)
        #rescaled and normalised x and y values (technically the xvals is not dependent on the data and maybe should be taken out of the loop for efficiency)
        xvals = [(x/xdivide)-zrange for x in range(zlength)]
        yvals = [float(zopen[x].strip())/(xdivide*10) for x in range(zstart,zstart+zlength)]
        
        #plt.figure("test1")
        #plt.scatter(xvals,yvals,s=5)
        #plt.savefig("test1" + str(i) + str(j) +".png")
        snippedindices = sorted(range(len(yvals)),key=lambda i: yvals[i])[-math.floor((topxpercent * 0.01) * len(yvals)):] 
        yvalssnipped = [yvals[k] for k in snippedindices]
        xvalssnipped = [xvals[k] for k in snippedindices]

        #plt.figure("test2")
        #plt.scatter(xvalssnipped,yvalssnipped,s=5)
        #plt.savefig("test2" + str(i) + str(j) +".png")
        # Create a random array to then randomise both yvals and xvals in the same way (this could be done by combining the arrays together, but I don't like being locked into that datastructure. I might do it if performance is bad down the line)
        sample = random.permutation(range(len(yvalssnipped)))
        # Splits up the x and y vals into 10 random subsets
        subsamplesx = [[xvalssnipped[x] for x in sample[math.floor((y * len(xvalssnipped)/subSampleAmtForGaussApprox)):math.floor(((y+1) * len(xvalssnipped))/subSampleAmtForGaussApprox)]] for y in range(subSampleAmtForGaussApprox)]
        subsamplesy = [[yvalssnipped[x] for x in sample[math.floor((y * len(yvalssnipped)/subSampleAmtForGaussApprox)):math.floor(((y+1) * len(xvalssnipped))/subSampleAmtForGaussApprox)]] for y in range(subSampleAmtForGaussApprox)]
        # Starting values for the gauss approx algorithm
        mean = sum([xvals[x]*yvals[x] for x in range(len(xvals))])/zlength
        sigma = sum([yvals[x]*(xvals[x]-mean)**2 for x in range(len(xvals))])/zlength
        # Loop over all subsamples
        for k in range(subSampleAmtForGaussApprox):
            #plt.figure("test")
            #plt.scatter(subsamplesx[k],subsamplesy[k],s=5)
            #plt.savefig("test" + str(i) + str(j) + str(k) + ".png")
            popt, pcov = curve_fit(gaus,subsamplesx[k],subsamplesy[k],p0=[0.2,mean,sigma])
            avgmean[i][j]+=popt[1]
            #print(popt)
            #print(popt[1])
            #print(minmean)
            #print(minmean[0])
            #print(minmean[0][0])
            if(popt[1] < minmean[i][j]):
                minmean[i][j] = popt[1]
            elif(popt[1] > maxmean[i][j]):
                maxmean[i][j] = popt[1]
        avgmean[i][j] = avgmean[i][j]/subSampleAmtForGaussApprox

print(len(avgmean))
for i in range(len(avgmean)):
    avgmean[i].insert(0,0)
print(len(avgmean))
print(len(avgmean[0]))
pushval.insert(0,"0")
plt.figure("rayplot")
meanerrormin = [[abs(avgmean[y][x] - minmean[y][x]) for x in range(len(pushval)-1)] for y in range(analysesteps)]
meanerrormax = [[abs(maxmean[y][x] - avgmean[y][x]) for x in range(len(pushval)-1)] for y in range(analysesteps)]
for i in range(len(meanerrormin)):
    meanerrormin[i].insert(0,0)
    meanerrormax[i].insert(0,0)
#np.transpose(meanerrormin)
#np.transpose(meanerrormax)
#avgmean = np.transpose(avgmean)

for i in range(analysesteps):
    plt.scatter([int(pushVal) for pushVal in pushval], avgmean[i], s=5)
    plt.errorbar([int(pushVal) for pushVal in pushval], avgmean[i], yerr=[meanerrormin[i],meanerrormax[i]],fmt="o")
plt.savefig("CC2D-" + str(size) + "-" + str(steps)+  "-rayplot.png")


nuval = [0. for x in range(analysesteps)]
for i in range(1,analysesteps):
    popt,pcov = curve_fit(line,[float('0.' + pushVal) for pushVal in pushval], avgmean[i])
    nuval[i] = math.log(2**analysesteps)/math.log(popt[0])
print(nuval)



        #    popt,pcov = curve_fit(gaus,xvals,yvals,p0=[1,mean,sigma])

#    plt.figure("z")
#    plt.scatter(xvals,yvals,s=5)
#    plt.figure("t")
#    plt.scatter(txvals,tvals,s=5)


#plt.figure("z")
#plt.savefig("CC2D-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-z.png")
#plt.figure("t")
#axes = plt.gca()
#axes.set_ylim([0,1000])
#plt.savefig("CC2D-" + str(size) + "-" + str(steps) + "-" + str(sym) + "-t.png")
