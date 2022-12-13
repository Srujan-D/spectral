import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os
pwd = os.getcwd()


input_txt = '/home/srujan/RISS/code/riss/src/btrapz/src/x_y_path.txt'
x = []
y = []

f = open(input_txt)

for line in f:
    line = line.strip('\n')
    line = line.split(' ')

    x.append(float(line[1]))
    y.append(float(line[0]))

f.close

# input_txt = '/home/ubuntu/scripts/piecewise_rectangle/plot_st.txt'
xr = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0,
	7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5]
yr = [0.0, 2.12, 3.0, 3.67, 4.24, 4.74, 5.2, 5.61, 6.0, 6.36, 6.71, 7.04, 7.35,
		  7.65, 7.94, 8.22, 8.49, 8.75, 9.0, 9.25, 9.49, 9.72, 9.95, 10.17, 10.39,
		  10.61, 10.82, 11.02, 11.22, 11.42]

# f = open(input_txt)

# for line in f:
#     line = line.strip('\n')
#     line = line.split(' ')

#     xr.append(float(line[0]))
#     yr.append(float(line[1]))

# f.close

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y, linewidth = 1.0, color = 'r', label = 'x y path')
ax1.legend(loc = 1)
ax1.plot(xr, yr, linewidth = 1.0,  color = 'b', label = 'reference x y path')

# ax1.plot(xr, yr, linewidth = 1.0, color = 'b', label = 'Rectangular Corridor')
ax1.legend(loc = 2)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
plt.title("X Y Path")

obs1_x = [4.0, 4.5, 7.5, 7]
obs1_y = [5, 7.5, 8.5, 6]

ax1.add_patch(patches.Polygon(xy=list(zip(obs1_x,obs1_y)), fill=True))

# obs2_x = [4.0, 6.0, 6.0, 4.0]
# obs2_y = [25, 45, 30, 10]

# ax1.add_patch(patches.Polygon(xy=list(zip(obs2_x,obs2_y)), fill=True))

plt.savefig('/home/srujan/RISS/code/riss/src/btrapz/src/figures/x_y_computed')

plt.show()
