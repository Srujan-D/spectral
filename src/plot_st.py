import math
from turtle import color
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
pwd = os.getcwd()

mscale = 10.0

input_txt = '/home/srujan_d/RISS/code/riss/src/btrapz/src/slt_3d.txt'
x = []
y = []

f = open(input_txt)

for line in f:
    line = line.strip('\n')
    line = line.split(' ')

    x.append(float(line[0]))
    y.append(float(line[1]))

f.close

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y, linewidth = 1.0, color = 'r', label = 'Trapezoidal Corridor')
ax1.legend(loc = 1)
ax1.set_xlabel("t / s")
ax1.set_ylabel("s / m")
plt.title("Piecewise Bezier Curve with Corridors")
obs1_y = [10, 16, 26, 20]
obs1_x = [10,12,12,10]

ax1.add_patch(patches.Polygon(xy=list(zip(obs1_x,obs1_y)), fill=True))

obs2_x = [15,17,17,15]
obs2_y = [12, 18, 28, 22]

ax1.add_patch(patches.Polygon(xy=list(zip(obs2_x,obs2_y)), fill=True))

corridor1_x = [9,9,10,10]
corridor1_y = [0,30,30,0]
ax1.add_patch(patches.Polygon(xy=list(zip(corridor1_x,corridor1_y)), fill=True, facecolor='orange', edgecolor='black'))

corridor1_x = [10,10,11,11]
corridor1_y = [0,10,13,0]
ax1.add_patch(patches.Polygon(xy=list(zip(corridor1_x,corridor1_y)), fill=True, facecolor='orange', edgecolor='black'))

corridor1_x = [11,11,12,12]
corridor1_y = [0,13,16,0]
ax1.add_patch(patches.Polygon(xy=list(zip(corridor1_x,corridor1_y)), fill=True, facecolor='orange', edgecolor='black'))

corridor1_x = [12,12,13,13]
corridor1_y = [0,30,30,0]
ax1.add_patch(patches.Polygon(xy=list(zip(corridor1_x,corridor1_y)), fill=True, facecolor='orange', edgecolor='black'))

corridor1_x = [13,13,14,14]
corridor1_y = [0,30,30,0]
ax1.add_patch(patches.Polygon(xy=list(zip(corridor1_x,corridor1_y)), fill=True, color='orange'))

plt.savefig('/home/srujan/RISS/code/riss/src/btrapz/src/figures/st_plot')

plt.show()