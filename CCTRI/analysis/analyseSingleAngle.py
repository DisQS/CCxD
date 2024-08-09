# to load matplotlib you need 
# module load GCC/11.2.0 OpenMPI/4.1.1
# module load matplotlib/3.5.2
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
from pathlib import Path


cur_path = Path.cwd()
new_path = cur_path.parent / "Data"

steps = 24
size = 6
spinangle = '0-0'
#spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000','2.000000','1.818000','1.667000','1.538000','1.429000','1.333000','1.250000','1.176400','1.111100','1.052600','1.000000']
#spinangles = ['20.000000','10.000000','6.666600','5.000000','3.333000','2.857000','2.500000','2.222000']
#spinangles = ['0.000000','0.050000','0.100000','0.150000','0.200000','0.250000','0.300000','0.350000','0.400000','0.450000']
thmax = []


print('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(spinangle) + '/' + str('1') + '/thdist.txt')
for i in range(1,24):
    th = open(new_path / ('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(spinangle)) / (str(i)) / 'thdist.txt', 'r').readlines()
    t = open(new_path /('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(spinangle)) / (str(i)) / 'tdist.txt', 'r').readlines()
    g = open(new_path /('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(spinangle)) / (str(i)) / 'gdist.txt', 'r').readlines()
    z = open(new_path /('CCTRI-'+str(size) + '-' + str(steps) + '-' + str(spinangle)) / (str(i)) / 'zdist.txt', 'r').readlines()
    
    dataxth = []
    datayth = []

    dataxt = []
    datayt = []

    dataxg = []
    datayg = []

    dataxz = []
    datayz = []

    for j in range(len(th)):
        dataxth.append(j/100)
        tmp = int(th[j].strip())
        datayth.append(tmp)
        

    for j in range(len(t)):
        dataxt.append(j/100)
        tmp = int(t[j].strip())
        datayt.append(tmp)

    for j in range(len(g)):
        dataxg.append(j/100)
        tmp = int(g[j].strip())
        datayg.append(tmp)

    for j in range(len(z)):
        dataxz.append((j/10)-25)
        tmp = int(z[j].strip())
        datayz.append(tmp)
    


    plt.figure("th")
    plt.scatter(dataxth,datayth,s=5)

    plt.figure("t")
    plt.scatter(dataxt,datayt,s=5)

    plt.figure("g")
    plt.scatter(dataxg,datayg,s=5)

    plt.figure("z")
    plt.scatter(dataxz,datayz,s=5)




plt.figure("th")
plt.savefig('thsa-'+str(size) + '-' + str(steps) + '-' + str(spinangle) + '.png')

plt.figure("t")
plt.savefig('tsa-'+str(size) + '-' + str(steps) + '-' + str(spinangle) + '.png')

plt.figure("g")
plt.savefig('gsa-'+str(size) + '-' + str(steps) + '-' + str(spinangle) + '.png')

plt.figure("z")
plt.savefig('zsa-'+str(size) + '-' + str(steps) + '-' + str(spinangle) + '.png')


