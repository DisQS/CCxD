# to load matplotlib you need 
# module load GCC/11.2.0 OpenMPI/4.1.1
# module load matplotlib/3.5.2
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker

steps = 40
size = 6
#spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000','2.000000','1.818000','1.667000','1.538000','1.429000','1.333000','1.250000','1.176400','1.111100','1.052600','1.000000']
#spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000']
#spinangles = ['0.000000','0.050000','0.100000','0.150000','0.200000','0.250000','0.300000','0.350000','0.400000','0.450000']
spinangles = ['0.000000','0.050000']
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

    ztotal = 0
    thtotal = 0
    gtotal = 0
    ttotal = 0

    for j in range(len(th)):
        dataxth.append(j/100)
        tmp = int(th[j].strip())
        datayth.append(tmp)
        thtotal+=tmp
        

    for j in range(len(t)):
        dataxt.append(j/100)
        tmp = int(t[j].strip())
        datayt.append(tmp)
        ttotal+=tmp

    for j in range(len(g)):
        dataxg.append(j/100)
        tmp = int(g[j].strip())
        datayg.append(tmp)
        gtotal+=tmp

    for j in range(len(z)):
        dataxz.append((j/20)-25)
        tmp = int(z[j].strip())
        datayz.append(tmp)
        ztotal+=tmp
    
    thmax.append(datayth.index(max(datayth)))
    znormal = [(20*z)/ztotal for z in datayz]
    tnormal = [(100*t)/ttotal for t in datayt]
    gnormal = [(100*g)/gtotal for g in datayg]
    thnormal = [(100*th)/thtotal for th in datayth]


    plt.figure("th")
    plt.scatter(dataxth,thnormal,s=5, label=str(float(i)))
    ax = plt.gca()
    ax.set_ylim([0, 3])
    plt.legend()

    plt.figure("t")
    plt.scatter(dataxt,tnormal,s=5)

    plt.figure("g")
    plt.scatter(dataxg,gnormal,s=5)

    plt.figure("z")
    plt.scatter(dataxz,znormal,s=5)


plt.figure("thmax")
plt.scatter([(float(i)) for i in spinangles],[i/100 for i in thmax],s=5)
plt.savefig("thmax.png")


plt.figure("th")
plt.savefig("th.png")

plt.figure("t")
plt.savefig("t.png")

plt.figure("g")
plt.savefig("g.png")

plt.figure("z")
plt.savefig("z.png")


