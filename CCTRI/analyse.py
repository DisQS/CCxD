import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker

steps = 10

for i in range(steps):
    th = open('outputth' + str(i) + '.txt', 'r').readlines()
    t = open('outputt' + str(i)  + '.txt', 'r').readlines()
    g = open('outputg' +str(i) + '.txt', 'r').readlines()
    z = open('outputz' + str(i) + '.txt', 'r').readlines()
    
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

    for j in range(len(z)):
        dataxz.append(j)
        datayz.append(int(z[j].strip()))

    plt.figure("th")
    plt.scatter(dataxth,datayth,s=5)

    plt.figure("t")
    plt.scatter(dataxt,datayt,s=5)

    plt.figure("g")
    plt.scatter(dataxg,datayg,s=5)

    plt.figure("z")
    plt.scatter(dataxz,datayz,s=5)

plt.figure("th")
plt.savefig("th.png")

plt.figure("t")
plt.savefig("t.png")

plt.figure("g")
plt.savefig("g.png")

plt.figure("z")
plt.savefig("z.png")


