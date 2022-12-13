import numpy as np

# scale an input array-like to a mininum and maximum number
# the input array must be of a floating point array
# if you have a non-floating point array, convert to floating using `astype('float')`
# this works with n-dimensional arrays
# it will mutate in place
# min and max can be integers
def scale_range (input, min, max):
    input += -(np.min(input))
    input /= np.max(input) / (max - min)
    input += min
    return input


input_txt = "/home/srujan_d/RISS/code/btrapz/src/s1_slt_3d_2.txt"
f = open(input_txt)
ego_y = []
ego_x = []
ego_dx = []
ego_dy = []

for line in f:
    line = line.strip("\n")
    line = line.split(" ")

    ego_y.append(float(line[1]))
    ego_x.append(float(line[2]))
    ego_dy.append(float(line[3]))
    ego_dx.append(float(line[4]))

f.close

# print(ego_x)
ego_x = np.array(ego_x)
np.interp(ego_x, (ego_x.min(), ego_x.max()), (0, 3.5))
ego_x.tolist()
# print(ego_x)


# print(scale_range(ego_x, 0.0, 3.5))

l = "1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.28 1.36 1.45 1.53 1.61 1.69 1.78 1.86 1.94 2.02 2.11 2.19 2.27 2.35 2.44 2.52 2.6 2.69 2.77 2.85 2.93 3.02 3.1 3.18 3.26 3.34 3.43 3.51 3.59 3.67 3.76 3.84 3.92 4.01 4.09 4.17 4.25 4.33 4.42 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5 4.5"

a = l.split(" ")
print(a)

# a.reverse()

# ego_x = np.array(a)
np.interp(ego_x, (ego_x.min(), ego_x.max()), (0, 3.5))
a = ego_x.tolist()

a.reverse()
print(a)

l = ''
for i in a:
    l += str(i) + ' '

print(l)

a = "0.0 6.15 0.0 6.6 0.0 7.05 0.0 7.5 0.0 7.95 0.0 8.4 0.0 8.85 0.0 9.3 0.0 9.75\
0.0 10.2 0.0 10.65 0.0 11.1 0.0 11.55 11.5 12.0 12.0 12.45 12.5 12.9 13.0 13.35 13.5 13.8 14.0 14.25 \
14.5 14.7 15.0 15.15 15.5 15.6 16.0 16.05 16.5 16.5 17.0 16.95 17.5 17.4 18.0 17.85 18.5 18.3 19.0 18.75 \
19.5 19.2 20.0 19.65 20.5 20.1 21.0 20.55 21.5 21.0 22.0 21.45 22.5 21.9 23.0 22.35 23.5 22.8 24.0 23.25 \
24.5 23.7 25.0 24.15 25.5 50.0 26.0 50.0 26.5 50.0 27.0 50.0 27.5 50.0 28.0 50.0 28.5 50.0 29.0 50.0 \
29.5 50.0 30.0 50.0 30.5 50.0 31.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 \
0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 0.0 50.0 \
0.0 50.0 0.0 50.0 "

b = a.split()
print(b)
print("\n\n\n")

a = ''
for i in range(71):
    a += "3.5 "

print(a)
print("\n\n\n")

b = [6*i/10 for i in range(71)]
a = ''
for i in range(71):
    a += str(round(b[i],2)) + " "
    if i%10 == 9:
        a+="\n"
print(a)

print("\n\n\n")

a = ''
for i in range(71):
    if i%10 == 9:
        a+="\n"
    if i < 11:
        a += "3.5 "
    elif i < 51:
        a += str(round(-3.75*i/40 + 50*3.75/40 - 0.25, 2)) + " "
    else:
        a += "-0.25 "

print(a)