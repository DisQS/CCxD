import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
import math
from pathlib import Path
import pylab as plb
from scipy.optimize import curve_fit
from numpy import random



steps = 32
size = 7
parallelamt = 0
analysesteps = 10
sym = 0 
zrange = 25
topxpercent = 10
subSampleAmtForGaussApprox = 10
thval = 0
phival = 0

errorfrombootstrap=0

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x,m,c):
    return (m * x) + c


cur_path = Path.cwd()
new_path = cur_path.parent /"Data"

foldername = 'CCTRI-st' + str(steps)+'-si'+str(size + parallelamt)+ ('-sym' if sym==1 else '') + '-th' + str(thval) + '-ph' + str(phival)
if not (Path(cur_path/foldername).is_dir()):
    Path.mkdir(cur_path / foldername)

pushval=["010","020","030","040","050","060","070"]
#pushval=["002","003","004","005","006","007","008"]
#pushval=["0010","0015","0020","0025","0030","0035","0040","0045","0050","0055","0060","0065","0070","0075","0080","0085","0090"]
z = [[None for y in range(analysesteps)] for x in range(len(pushval))]
avgmean = [[0.0 for y in range(len(pushval))] for x in range(analysesteps)]
if not errorfrombootstrap:
   minmean = [[50.0 for y in range(len(pushval))] for x in range(analysesteps)] 
   maxmean = [[-50.0 for y in range(len(pushval))] for x in range(analysesteps)]
errormean = [[0.0 for y in range(len(pushval))] for x in range(analysesteps)]



print("Importing and Analysing data")
for i in range(0,analysesteps):
    for j in range(0,len(pushval)):
        zopen = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + "0-0") / pushval[j] / ("dists" + str(i+1) + ".txt" )).readlines()

        print(  (str(math.floor(((i) * len(pushval) + (j+1))/(len(pushval) * analysesteps)*100))) + '%', end='\r' )
        #slog through the data and find out where the z values actually start

        thlength = int(zopen[0])
        tlength = int(zopen[1])
        glength = int(zopen[2])
        zlength = int(zopen[3])
        zstart = 4+int(thlength)+int(tlength)+int(glength)
        
        #constant to rescale and normalise data
        zxdivide = zlength/(2*zrange)
        zydivide = (10**size)/zxdivide
        #rescaled and normalised x and y values (technically the xvals is not dependent on the data and maybe should be taken out of the loop for efficiency)
        xvals = [(x/zxdivide)-zrange for x in range(zlength)]
        yvals = [float(zopen[x].strip())/(zydivide) for x in range(zstart,zstart+zlength)]
        
        
        snippedindices = sorted(range(len(yvals)),key=lambda i: yvals[i])[-math.floor((topxpercent * 0.01) * len(yvals)):] 
        yvalssnipped = [yvals[k] for k in snippedindices]
        xvalssnipped = [xvals[k] for k in snippedindices]

        if not errorfrombootstrap:
            # Create a random array to then randomise both yvals and xvals in the same way (this could be done by combining the arrays together, but I don't like being locked into that datastructure. I might do it if performance is bad down the line)
            sample = random.permutation(range(len(yvalssnipped)))
            # Splits up the x and y vals into 10 random subsets
            subsamplesx = [[xvalssnipped[x] for x in sample[math.floor((y * len(xvalssnipped)/subSampleAmtForGaussApprox)):math.floor(((y+1) * len(xvalssnipped))/subSampleAmtForGaussApprox)]] for y in range(subSampleAmtForGaussApprox)]
            subsamplesy = [[yvalssnipped[x] for x in sample[math.floor((y * len(yvalssnipped)/subSampleAmtForGaussApprox)):math.floor(((y+1) * len(xvalssnipped))/subSampleAmtForGaussApprox)]] for y in range(subSampleAmtForGaussApprox)]
        # Starting values for the gauss approx algorithm
        mean = sum([xvals[x]*yvals[x] for x in range(len(xvals))])/zlength
        sigma = sum([yvals[x]*(xvals[x]-mean)**2 for x in range(len(xvals))])/zlength
        # Loop over all subsamples

        if not errorfrombootstrap:
            #Error bars with bootstrapping
            for k in range(subSampleAmtForGaussApprox):
            
                popt, pcov = curve_fit(gaus,subsamplesx[k],subsamplesy[k],p0=[0.2,mean,sigma])
                avgmean[i][j]+=popt[1]
                if(popt[1] < minmean[i][j]):
                    minmean[i][j] = popt[1]
                elif(popt[1] > maxmean[i][j]):
                    maxmean[i][j] = popt[1]
            avgmean[i][j] = avgmean[i][j]/subSampleAmtForGaussApprox

        else:
            #errorbars with curve_fit's own error
            popt, pcov = curve_fit(gaus,xvalssnipped, yvalssnipped, p0=[0.2,mean,sigma])
            perr = np.sqrt((np.diag(pcov)))
            avgmean[i][j] = popt[1]
            errormean[i][j] = perr[1]


print("Done!")

for i in range(len(avgmean)):
    avgmean[i].insert(0,0)

pushval.insert(0,"0")
rayfig, rayax = plt.subplots()
rayax.set_xlabel(r'$z^0_\text{max}$')
rayax.set_ylabel(r'$z^n_\text{max}$')
rayax.set_title('Rayplot for 10^' + str(size) + ' samples, ' + str(steps) + ' steps')
rayax.set_axisbelow(True)
rayax.grid(True)
rayax.set_xlim(0.0,1.1 * max([int(pushVal)/1000 for pushVal in pushval]))
rayax.set_ylim(0.0,1.2 * avgmean[-1][-1])


if not errorfrombootstrap:
    meanerrormin = [[abs(avgmean[y][x] - minmean[y][x]) for x in range(len(pushval)-1)] for y in range(analysesteps)]
    meanerrormax = [[abs(maxmean[y][x] - avgmean[y][x]) for x in range(len(pushval)-1)] for y in range(analysesteps)]
for i in range(len(errormean)):
    errormean[i].insert(0,0)
    if not errorfrombootstrap:
        meanerrormin[i].insert(0,0)
        meanerrormax[i].insert(0,0)

for i in range(analysesteps):
    rayax.scatter([int(pushVal)/1000 for pushVal in pushval], avgmean[i], s=5)
    rayax.errorbar([int(pushVal)/1000 for pushVal in pushval], avgmean[i], yerr=[errormean[i],errormean[i]],fmt="o")
    popt,pcov = curve_fit(line,[int(pushVal)/1000 for pushVal in pushval], avgmean[i])
    linevals = np.arange(0.0,1.1 * int(pushval[-1])/1000 , 0.0001)
    lineyvals = popt[1] + popt[0] * linevals
    rayax.plot(linevals,lineyvals)
rayfig.savefig(foldername + "/CCTRI-" + str(size) + "-" + str(steps)+"-" + str(thval) + "-" + str(phival)+  "-rayplot.png")


nuval = [0. for x in range(analysesteps)]
nuerror = [0. for x in range(analysesteps)]
for i in range(1,analysesteps+1):
    popt,pcov = curve_fit(line,[int(pushVal)/1000 for pushVal in pushval], avgmean[i-1])
    print(math.log(2**i))
    print(math.log(popt[0]))
    nuval[i-1] = math.log(2**i)/math.log(popt[0])
    nuerror[i-1] = np.sqrt(np.diag(pcov))[0]
nuplotxval = [2**n for n in range(1,analysesteps+1)]

nufig, nuax = plt.subplots()
nuax.set_xlabel('$2^n$')
nuax.set_ylabel('$ν$')
nuax.set_title("critical exponent graph")
nuax.set_axisbelow(True)
nuax.grid(True)

nuax.scatter(nuplotxval, nuval, s=5)
nuax.errorbar(nuplotxval,nuval,yerr=[nuerror,nuerror],fmt="o")
nufig.savefig(foldername + "/CCTRI-" +str(size) + "-" + str(steps) + "-" + str(thval) + "-" + str(phival) + "-nuplot.png")
print(nuval)
print(nuerror)
