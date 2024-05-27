import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as maticker
import math
from pathlib import Path

steps = 40
size = 6
sym = 0

cur_path = Path.cwd()
new_path = cur_path.parent / ("TRIRG-"  + str(size) + "-" + str(steps) + "-" + str(sym)) / "Data"



thmin = 0
thmax = 40
thstep = 2

psimin = 0
psimax = 40
psistep = 2

vals = [[0 for x in range(math.ceil((thmax-thmin)/thstep))] for y in range(math.ceil((psimax-psimin)/psistep))]

for th in range(thmin,thmax,thstep):

    for psi in range(psimin,psimax,psistep):

        i = math.floor((th-thmin)/thstep)
        j = math.floor((psi-psimin)/psistep)
        try:
            z = open(new_path / ("CCTRI-" + str(size) + "-" + str(steps) + "-" + str(th) + "-" + str(psi)) / str(steps) / "zdist.txt" ).readlines()
            break
        except FileNotFoundError:
            z = [0]
        
        
        for val in range(len(z)):
            vals[i][j]+=int(z[val].strip())
        
    
print(vals)

plt.figure("phasemap")
plt.imshow(vals,cmap="hot",interpolation="nearest")
plt.savefig("phasemap.png")
