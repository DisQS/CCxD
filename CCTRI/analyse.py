# to load matplotlib you need 
# module load GCC/11.2.0 OpenMPI/4.1.1
# module load matplotlib/3.5.2
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker

steps = 20
size = 6
#spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000','2.000000','1.818000','1.667000','1.538000','1.429000','1.333000','1.250000','1.176400','1.111100','1.052600','1.000000']
spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000']
thmax = []


for i in spinangles:
    th = open('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(i) + '/' + str(steps-1) + '/thdist.txt', 'r').readlines()
    t = open('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(i) + '/' + str(steps-1) + '/tdist.txt', 'r').readlines()
    g = open('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(i) + '/' + str(steps-1) + '/gdist.txt', 'r').readlines()
    z = open('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(i) + '/' + str(steps-1) + '/zdist.txt', 'r').readlines()
    
    dataxth = []
    datayth = []

    dataxt = []
    datayt = []

    dataxg = []
    datayg = []

    dataxz = []
    datayz = []

    for j in range(len(th)):
        dataxth.append(j)
        datayth.append(int(th[j].strip()))

    for j in range(len(t)):
        dataxt.append(j)
        datayt.append(int(t[j].strip()))

    for j in range(len(g)):
        dataxg.append(j)
        datayg.append(int(g[j].strip()))

    for j in range(1,len(z)):
        dataxz.append(j)
        datayz.append(int(z[j].strip()))
    
    thmax.append(datayth.index(max(datayth)))

    plt.figure("th")
    plt.scatter(dataxth,datayth,s=5, label=str(2/float(i)))
    plt.legend()

    plt.figure("t")
    plt.scatter(dataxt,datayt,s=5)

    plt.figure("g")
    plt.scatter(dataxg,datayg,s=5)

    plt.figure("z")
    plt.scatter(dataxz,datayz,s=5)


plt.figure("thmax")
plt.scatter([(2/float(i)) for i in spinangles],[i/100 for i in thmax],s=5)
plt.savefig("thmax.png")


plt.figure("th")
plt.savefig("th.png")

plt.figure("t")
plt.savefig("t.png")

plt.figure("g")
plt.savefig("g.png")

plt.figure("z")
plt.savefig("z.png")


